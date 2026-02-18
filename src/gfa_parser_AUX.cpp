#include "../include/gfa_parser_AUX.hpp"

#include <sstream>    // for std::stringstream, std::ostringstream


const uint8_t* GfaAuxParser::skip_tag_data(const uint8_t* cur, const uint8_t* end) {
    if (cur >= end) return nullptr;
    char type = static_cast<char>(*cur++);
    switch (type) {
        case 'A':  return cur + 1;
        case 'i':
        case 'f':  return cur + 4;
        case 'Z':  while (cur < end && *cur) ++cur; return (cur < end)? cur + 1 : nullptr;
        case 'B': {
            if (cur + 5 > end) return nullptr;
            char subtype = *cur++;
            (void)subtype;
            int32_t n = *(reinterpret_cast<const int32_t*>(cur)); cur += 4;
            size_t elem_sz = 0;
            switch (subtype) {
                case 'c': case 'C': elem_sz = 1; break;
                case 's': case 'S': elem_sz = 2; break;
                case 'i': case 'I': case 'f': elem_sz = 4; break;
                default: return nullptr;
            }
            size_t bytes = static_cast<size_t>(n) * elem_sz;
            if (end - cur < static_cast<ptrdiff_t>(bytes)) return nullptr;
            return cur + bytes;
        }
        default: return nullptr;
    }
}

std::vector<uint8_t> GfaAuxParser::parse_aux_string(const std::string& s_aux) {
    std::vector<uint8_t> out;
    if (s_aux.empty()) { return out; }

    std::string str = s_aux;
    if (!str.empty() && str[0] == '\t') str.erase(0,1);

    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, '\t')) {
        if (token.size() < 5 || token[2] != ':' || token[4] != ':') continue;
        char k1 = token[0], k2 = token[1], typ = token[3];
        std::string val = token.substr(5);

        out.push_back(static_cast<uint8_t>(k1));
        out.push_back(static_cast<uint8_t>(k2));
        out.push_back(static_cast<uint8_t>(typ));

        auto push32 = [&](uint32_t v){ const uint8_t* p=reinterpret_cast<uint8_t*>(&v); out.insert(out.end(),p,p+4); };
        auto pushfloat = [&](float f){ const uint8_t* p=reinterpret_cast<uint8_t*>(&f); out.insert(out.end(),p,p+4); };

        switch (typ) {
            case 'A': out.push_back(val.empty()? 0 : static_cast<uint8_t>(val[0])); break;
            case 'i': push32(static_cast<int32_t>(std::stol(val))); break;
            case 'f': pushfloat(std::stof(val)); break;
            case 'Z': out.insert(out.end(), val.begin(), val.end()); out.push_back(0); break;
            case 'B':  /* store as Z for simplicity */ 
                out.back() = static_cast<uint8_t>('Z');
                out.insert(out.end(), val.begin(), val.end()); out.push_back(0); break;
            default:  /* unknown -> store as Z */ 
                out.back() = static_cast<uint8_t>('Z');
                out.insert(out.end(), val.begin(), val.end()); out.push_back(0); break;
        }
    }
    return out;
}

const uint8_t* GfaAuxParser::get_aux_data(const std::vector<uint8_t>& aux_data, const char tag[2])
{
    const uint8_t* cur = aux_data.data();
    const uint8_t* end = aux_data.data() + aux_data.size();
    uint16_t want = static_cast<uint16_t>(tag[0] << 8 | tag[1]);
    while (cur && cur + 3 <= end) {
        uint16_t key = static_cast<uint16_t>(cur[0] << 8 | cur[1]);
        const uint8_t* type_ptr = cur + 2;
        cur = type_ptr;
        if (key == want) return type_ptr;
        cur = skip_tag_data(cur, end);
    }
    return nullptr;
}

void GfaAuxParser::delete_aux_tag(std::vector<uint8_t>& aux_data, const uint8_t* type_ptr) {
    if (!type_ptr) return;
    const uint8_t* tag_begin = type_ptr - 2;
    const uint8_t* tag_end   = skip_tag_data(type_ptr, aux_data.data() + aux_data.size());
    if (!tag_end) return;
    aux_data.erase(
        aux_data.begin() + (tag_begin - aux_data.data()),
        aux_data.begin() + (tag_end   - aux_data.data())
    );
}

void GfaAuxParser::append_int_tag(std::vector<uint8_t>& aux_data, const char tag[2], int32_t val) {
    // format: 2 bytes key + 1 byte type 'i' + 4 bytes value (little-endian)
    aux_data.push_back(static_cast<uint8_t>(tag[0]));
    aux_data.push_back(static_cast<uint8_t>(tag[1]));
    aux_data.push_back(static_cast<uint8_t>('i'));
    const uint8_t* p = reinterpret_cast<const uint8_t*>(&val);
    aux_data.insert(aux_data.end(), p, p+4);
}

void GfaAuxParser::append_str_tag(std::vector<uint8_t>& aux_data, const char tag[2], const std::string& s) {
    // format: 2 bytes key + 1 byte type 'Z' + string value + null terminator
    aux_data.push_back(static_cast<uint8_t>(tag[0]));
    aux_data.push_back(static_cast<uint8_t>(tag[1]));
    aux_data.push_back(static_cast<uint8_t>('Z'));
    aux_data.insert(aux_data.end(), s.begin(), s.end());
    aux_data.push_back(0);
}

std::string GfaAuxParser::write_aux_as_text(const std::vector<uint8_t>& aux_data) {
    std::ostringstream oss;

    const uint8_t* cur = aux_data.data();
    const uint8_t* end = aux_data.data() + aux_data.size();
    while (cur && cur + 3 <= end) {
        uint8_t k1 = cur[0], k2 = cur[1];
        const uint8_t* type_ptr = cur + 2;
        cur = type_ptr;

        if (cur >= end) break;
        char typ = static_cast<char>(*cur++);
        oss << '\t' << char(k1) << char(k2) << ':' << typ << ':';

        if (typ == 'A') {
            if (cur < end) oss << char(*cur++);
        } else if (typ == 'i') {
            if (cur + 4 <= end) { oss << read_le_i32(cur); cur += 4; }
        } else if (typ == 'f') {
            if (cur + 4 <= end) { oss << read_le_f32(cur); cur += 4; }
        } else if (typ == 'Z') {
            const uint8_t* s = cur;
            while (s < end && *s) ++s;
            oss.write(reinterpret_cast<const char*>(cur), s - cur);
            cur = (s < end ? s + 1 : s);
        } else if (typ == 'B') {
            if (cur + 5 <= end) {
                char subtype = char(*cur++);
                int32_t n = read_le_i32(cur); cur += 4;
                size_t elem_sz = (subtype=='c'||subtype=='C')?1:(subtype=='s'||subtype=='S')?2:4;
                size_t bytes = size_t(n) * elem_sz;
                oss << subtype << ',' << n << ',';

                const uint8_t* p = cur;
                for (int i = 0; i < n; ++i) {
                    if (i) oss << ',';
                    if (elem_sz == 1) oss << int(*p), ++p;
                    else if (elem_sz == 2) {
                        std::int16_t v16 = 0;
                        std::copy_n(p, 2, reinterpret_cast<std::uint8_t*>(&v16));
                        oss << v16;
                        p += 2;
                    }
                    else {
                        std::int32_t v32 = 0;
                        std::copy_n(p, 4, reinterpret_cast<std::uint8_t*>(&v32));
                        oss << v32;
                        p += 4;
                    }
                }
                cur += bytes;
            } else {
                // Tolerant: treat as empty Z
            }
        } else {
            // Unknown type: skip without crashing
            const uint8_t* nxt = GfaAuxParser::skip_tag_data(type_ptr, end);
            cur = nxt ? nxt : end;
        }
    }
    return oss.str();
}
#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>

// AUX tag (such as LN, SR, etc.) parser
struct GfaAux {
    std::vector<uint8_t> aux_data;  // [ key1 (1B) | key2 (1B) | type (1B) | value (N B) ] : LN:i:6 -> [L | N | i | 6 (4B)]
};

// AUX tag parser
class GfaAuxParser {
public:

    /**
     * @brief skip over value bytes
     *
     * @param cur: pointer to the type byte
     * @param end: pointer to the end of the aux_data vector
     * @return: pointer to the byte after the value, or nullptr if invalid
     */
    static const uint8_t* skip_tag_data(const uint8_t* cur, const uint8_t* end);

    /**
     * @brief parse auxiliary string into raw bytes
     *
     * @param s_aux: auxiliary string to parse
     * @return: vector of raw bytes
     */
    static std::vector<uint8_t> parse_aux_string(const std::string& s_aux);

    /**
     * @brief retrieve pointer to type byte
     *
     * @param aux_data: vector of raw bytes
     * @param tag: 2-byte tag to search for
     * @return: pointer to the type byte, or nullptr if not found
     */
    static const uint8_t* get_aux_data(const std::vector<uint8_t>& aux_data, const char tag[2]);

    /**
     * @brief delete auxiliary tag from raw bytes
     *
     * @param aux_data: vector to modify
     * @param type_ptr: pointer to the type byte of the tag to delete
     */
    static void delete_aux_tag(std::vector<uint8_t>& aux_data, const uint8_t* ptr_to_type_byte);

    /**
     * @brief append integer tag
     * 
     * @param aux_data: vector to modify
     * @param tag: 2-byte tag to append
     * @param val: integer value to append
     */
    static void append_int_tag(std::vector<uint8_t>& aux_data, const char tag[2], int32_t val);

    /**
     * @brief append string tag
     * 
     * @param aux_data: vector to modify
     * @param tag: 2-byte tag to append
     * @param s: string value to append
     */
    static void append_str_tag(std::vector<uint8_t>& aux_data, const char tag[2], const std::string& s);

    /**
     * @brief write auxiliary data as text
     *
     * @param aux_data: vector of raw bytes
     * @return: formatted string representation of the auxiliary data
     */    
    static std::string write_aux_as_text(const std::vector<uint8_t>& aux_data);

private:

    static inline int32_t read_le_i32(const uint8_t* p) {
        int32_t v;
        std::copy_n(p, 4, reinterpret_cast<std::uint8_t*>(&v));
        return v;
    }
    static inline float read_le_f32(const uint8_t* p) {
        float f;
        std::copy_n(p, 4, reinterpret_cast<std::uint8_t*>(&f));
        return f;
    }

};
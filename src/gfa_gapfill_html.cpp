#include "gfa_gapfill_html.hpp"
#include "gfa_gapfill.hpp"

#include "logger.hpp"
#include "save.hpp"

#include <iomanip>
#include <sstream>

namespace {

std::string json_string(const std::string& value) {
    std::string out;
    out.reserve(value.size() + 2);
    out.push_back('"');
    for (unsigned char c : value) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    constexpr char hex[] = "0123456789abcdef";
                    out += "\\u00";
                    out.push_back(hex[c >> 4]);
                    out.push_back(hex[c & 15]);
                } else {
                    out.push_back(static_cast<char>(c));
                }
        }
    }
    out.push_back('"');
    return out;
}

}  // namespace

void GfaGapfill::write_html_(const std::vector<Candidate>& candidates) const {
    if (params_.html_file.empty()) return;

    std::vector<gapfill_html::Contig> contigs;
    contigs.reserve(fragments_.size());
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible) continue;
        contigs.push_back({
            i, fragment.component, fragment.hap, fragment.length,
            fragment.layout_start, fragment.reverse,
            paths_[fragment.path_id].name, fragment.sample,
            nodes_[Vertex::get_segment_id(fragment.vertices.front())].name +
                (Vertex::get_is_reverse(fragment.vertices.front()) ? "-" : "+"),
            nodes_[Vertex::get_segment_id(fragment.vertices.back())].name +
                (Vertex::get_is_reverse(fragment.vertices.back()) ? "-" : "+")
        });
    }

    std::vector<gapfill_html::Connection> connections;
    connections.reserve(candidates.size());
    for (const Candidate& candidate : candidates) {
        gapfill_html::Connection connection;
        connection.component = fragments_[candidate.left].component;
        connection.left = candidate.left;
        connection.right = candidate.right;
        connection.bridge = candidate.bridge;
        connection.phase_left = candidate.left_boundary.phase.similarity;
        connection.phase_right = candidate.right_boundary.phase.similarity;
        connection.phase_left_target_bp = candidate.left_boundary.phase.target_bp;
        connection.phase_left_bridge_bp = candidate.left_boundary.phase.bridge_bp;
        connection.phase_right_target_bp = candidate.right_boundary.phase.target_bp;
        connection.phase_right_bridge_bp = candidate.right_boundary.phase.bridge_bp;
        connection.boundary_left = candidate.left_boundary.minimizer;
        connection.boundary_right = candidate.right_boundary.minimizer;
        connection.identity_left = candidate.left_boundary.identity;
        connection.coverage_left = candidate.left_boundary.coverage;
        connection.identity_right = candidate.right_boundary.identity;
        connection.coverage_right = candidate.right_boundary.coverage;
        connection.sample_score = candidate.sample_score;
        connection.phase_score = candidate.phase_score;
        connection.alignment_score = candidate.alignment_score;
        connection.confidence = candidate.confidence;
        connection.homolog_span = candidate.homolog_span;
        connection.sample_support = candidate.sample_support;
        connection.spanning_samples = candidate.spanning_samples;

        const bool selected = candidate.status == Candidate::KEPT;
        connection.left_cut_bp = selected ? candidate.left_boundary.target_cut_bp :
            fragments_[candidate.left].path_bp[candidate.left_boundary.target_pos + 1];
        connection.right_cut_bp = selected ? candidate.right_boundary.target_cut_bp :
            fragments_[candidate.right].path_bp[candidate.right_boundary.target_pos];
        connection.bridge_left_bp = selected ? candidate.left_boundary.bridge_cut_bp :
            fragments_[candidate.bridge].path_bp[candidate.left_boundary.bridge_pos + 1];
        connection.bridge_right_bp = selected ? candidate.right_boundary.bridge_cut_bp :
            fragments_[candidate.bridge].path_bp[candidate.right_boundary.bridge_pos];

        switch (candidate.status) {
            case Candidate::KEPT: connection.status = "keep"; break;
            case Candidate::LOW_CONFIDENCE: connection.status = "drop:low-confidence"; break;
            case Candidate::USED_END: connection.status = "drop:used-end"; break;
            case Candidate::PHASE_CONFLICT: connection.status = "drop:phase-conflict"; break;
            case Candidate::COORDINATE_CONFLICT: connection.status = "drop:coordinate-conflict"; break;
            case Candidate::CYCLE: connection.status = "drop:cycle"; break;
            default: connection.status = "pending"; break;
        }
        connections.push_back(std::move(connection));
    }
    gapfill_html::save(params_.html_file, contigs, connections);
}

void gapfill_html::save(
    const std::string& file,
    const std::vector<Contig>& contigs,
    const std::vector<Connection>& connections
) {
    log_stream() << "Writing interactive gap-fill candidate graph ...\n";

    SAVE out(file);
    out.save(R"HTML(<!doctype html>
<html><head><meta charset="utf-8"><title>liftasm gap-fill candidates</title>
<style>
*{box-sizing:border-box}body{margin:0;font:14px system-ui,sans-serif;color:#182338;background:#f6f8fc}
header{height:54px;padding:10px 18px;background:linear-gradient(105deg,#122039,#29466d);color:white;font-size:18px;font-weight:650;letter-spacing:.2px;box-shadow:0 2px 8px #10192a33}
header small{margin-left:10px;color:#c9d7ea;font-size:12px;font-weight:450}
main{display:grid;grid-template-columns:400px 1fr;height:calc(100vh - 54px)}
aside{display:flex;min-height:0;flex-direction:column;border-right:1px solid #d8dee9;background:#fff;overflow:hidden;box-shadow:2px 0 8px #1720330a;z-index:2}.controls{padding:14px 16px 7px;border-bottom:1px solid #e4e9f0;background:#fbfcfe}
label{display:block;margin:0 0 5px;color:#536174;font-size:12px;font-weight:700;text-transform:uppercase}
select,input[type=search]{width:100%;padding:7px;margin-bottom:12px;border:1px solid #b8c1cf;border-radius:5px;background:white}input[type=range]{width:100%;margin:2px 0 14px}
#selection{padding:12px 16px;border-bottom:1px solid #e1e6ee;background:#fff;flex:none}#selection.focused{background:linear-gradient(135deg,#fff8f4,#fff)}#selection.focused label{color:#c94622}#detail{max-height:34vh;overflow:auto;white-space:pre-wrap;padding:11px 12px;border:1px solid #dfe5ed;border-radius:8px;background:#f4f6f9;line-height:1.48;overflow-wrap:anywhere}#selection.focused #detail{border-color:#f0a48e;background:#fff;box-shadow:0 4px 14px #c9462215}
.coordinates{margin-top:9px}.evidence{width:100%;margin:10px 0 8px;border-collapse:collapse;font-size:11px;white-space:nowrap}.evidence th,.evidence td{padding:5px 4px;border-bottom:1px solid #e7eaf0;text-align:right}.evidence th:first-child{text-align:left}.evidence thead th{color:#607087;font-weight:700}.support{color:#536174;font-size:12px}
#candidatePanel{min-height:0;flex:1;overflow:auto;padding:11px 16px 14px}#summary{margin:0 0 10px;color:#536174}.tools{display:flex;gap:7px;margin:-4px 0 7px}.tool{padding:6px 11px;border:1px solid #bcc6d4;border-radius:6px;background:#fff;color:#34435a;cursor:pointer}.tool:hover{background:#edf4ff}
.candidate{width:100%;padding:10px 11px;margin:0 0 8px;text-align:left;border:1px solid #dde3ec;border-left:4px solid #aab4c3;border-radius:8px;background:white;cursor:pointer;transition:background .14s,border-color .14s,box-shadow .14s,transform .14s}.candidate.keep{border-left-color:#159447}.candidate:hover{background:#edf4ff;transform:translateY(-1px);box-shadow:0 4px 10px #1b35551a}.candidate.active{border-color:#e4572e;border-left:7px solid #e4572e;background:linear-gradient(110deg,#fff0e9,#fff);box-shadow:0 0 0 2px #e4572e26,0 7px 18px #b83f1e24;transform:none}.candidate.active .route{color:#b53618}.route{display:block;font-weight:700;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.via,.metrics{display:block;margin-top:3px;color:#657286;font-size:11px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
#view{overflow:auto;position:relative;background:#fff;scroll-behavior:smooth}svg{min-height:100%;background:#fff}.lane{fill:#f8faff}.lane.alt{fill:#fff}.row{stroke:#e2e7ef;stroke-width:1}.axis{font-size:10px;fill:#77859a}.sample{font-size:12px;font-weight:700;fill:#34435a}.edge{fill:none;stroke:#9aa5b5;stroke-width:1.5;opacity:.20;cursor:pointer;transition:opacity .18s,stroke-width .18s}.edge.keep{stroke:#159447;stroke-width:2.5;opacity:.82}.edge.related{opacity:.42}.edge.active{stroke:#e4572e;stroke-width:4;opacity:1}.edge.dim{opacity:.025}.bridge-span{stroke-width:7;stroke-linecap:round;opacity:.32}.bridge-span.keep{opacity:.9}.anchor{stroke:white;stroke-width:1.5;fill:#7d8999;opacity:.7}.anchor.keep{fill:#159447;opacity:1}.anchor.active{fill:#e4572e}
.node{transition:opacity .18s}.node rect{fill:#edf2f8;fill-opacity:.82;stroke:#8492a6;stroke-width:1}.node.h1 rect{fill:#d9eaff;stroke:#3274c5}.node.h2 rect{fill:#ffe2e5;stroke:#d25763}.node.active rect{fill-opacity:1;stroke:#e4572e;stroke-width:3;filter:drop-shadow(0 3px 5px #e4572e55)}.node text{font-size:10px;pointer-events:none;fill:#1d2a3d}
.legend{display:flex;gap:12px;margin:2px 0 12px;color:#536174}.dot:before{content:'';display:inline-block;width:12px;height:3px;margin-right:5px;vertical-align:middle;background:#9aa5b5}.dot.keep:before{background:#159447}
#tooltip{position:fixed;display:none;max-width:420px;padding:7px 9px;border-radius:5px;background:#17243a;color:white;font-size:12px;line-height:1.4;white-space:pre-line;pointer-events:none;z-index:5;box-shadow:0 3px 12px #10182855}
</style></head><body><header>liftasm gap-fill <small>contig → bridge interval → contig</small></header><main><aside><div class="controls">
<label for="component">Component</label><select id="component"></select>
<label for="status">Connections</label><select id="status"><option value="all">All candidates</option><option value="keep">Kept only</option><option value="drop">Rejected only</option></select>
<label for="query">Find candidate</label><input id="query" type="search" placeholder="contig or sample">
<label for="zoom">Horizontal zoom: <span id="zoomValue">1x</span></label><input id="zoom" type="range" min="1" max="20" value="1" step="1">
<div class="tools"><button class="tool" id="fit" type="button">Fit view</button><button class="tool" id="clear" type="button">Clear selection</button></div>
<div class="legend"><span class="dot keep">kept</span><span class="dot">rejected</span></div></div>
<section id="selection"><label>Selected candidate</label><div id="detail">Click a connection or candidate.</div></section>
<section id="candidatePanel"><div id="summary"></div><div id="list"></div></section>
</aside><div id="view"><svg id="graph"></svg></div></main><div id="tooltip"></div><script>
const N=)HTML");

    out.save("[");
    for (size_t i = 0; i < contigs.size(); ++i) {
        const Contig& contig = contigs[i];
        if (i) out.save(",");
        std::ostringstream line;
        line << "{i:" << contig.id << ",c:" << contig.component
             << ",h:" << static_cast<uint32_t>(contig.hap)
             << ",b:" << contig.bp << ",x:" << contig.start
             << ",rv:" << (contig.reverse ? 1 : 0)
             << ",n:" << json_string(contig.name)
             << ",s:" << json_string(contig.sample)
             << ",a:" << json_string(contig.first_vertex)
             << ",z:" << json_string(contig.last_vertex) << "}";
        out.save(line.str());
    }
    out.save("];\nconst E=[");
    for (size_t i = 0; i < connections.size(); ++i) {
        const Connection& connection = connections[i];
        if (i) out.save(",");
        std::ostringstream line;
        line << std::setprecision(8)
             << "{i:" << i << ",c:" << connection.component
             << ",l:" << connection.left << ",r:" << connection.right
             << ",b:" << connection.bridge
             << ",pl:" << connection.phase_left << ",prh:" << connection.phase_right
             << ",plt:" << connection.phase_left_target_bp
             << ",plb:" << connection.phase_left_bridge_bp
             << ",prt:" << connection.phase_right_target_bp
             << ",prb:" << connection.phase_right_bridge_bp
             << ",ml:" << connection.boundary_left
             << ",mr:" << connection.boundary_right
             << ",il:" << connection.identity_left
             << ",cl:" << connection.coverage_left
             << ",ir:" << connection.identity_right
             << ",cr:" << connection.coverage_right
             << ",sscore:" << connection.sample_score
             << ",qscore:" << connection.phase_score
             << ",ascore:" << connection.alignment_score
             << ",cf:" << connection.confidence
             << ",ss:" << connection.sample_support
             << ",sp:" << connection.spanning_samples
             << ",hs:" << (connection.homolog_span ? 1 : 0)
             << ",lp:" << connection.left_cut_bp
             << ",rp:" << connection.right_cut_bp
             << ",bl:" << connection.bridge_left_bp
             << ",br:" << connection.bridge_right_bp
             << ",st:" << json_string(connection.status) << "}";
        out.save(line.str());
    }
    out.save(R"HTML(];
const $=function(id){return document.getElementById(id)},svg=$('graph'),view=$('view'),comp=$('component'),filter=$('status'),query=$('query'),zoom=$('zoom'),zoomValue=$('zoomValue'),list=$('list'),detail=$('detail'),selection=$('selection'),summary=$('summary'),tooltip=$('tooltip');
const ns='http://www.w3.org/2000/svg',byId=new Map(N.map(function(n){return[n.i,n]}));
function S(tag,attrs,text){const x=document.createElementNS(ns,tag);Object.keys(attrs||{}).forEach(function(k){x.setAttribute(k,attrs[k])});if(text!==undefined)x.textContent=text;return x}
function visible(e){const status=filter.value==='all'||(filter.value==='keep'&&e.st==='keep')||(filter.value==='drop'&&e.st!=='keep');if(!status)return false;const q=query.value.trim().toLowerCase();if(!q)return true;const a=byId.get(e.l),b=byId.get(e.r),c=byId.get(e.b);return[a.n,b.n,c.n,a.s,b.s,c.s,e.st].join(' ').toLowerCase().includes(q)}
function fmt(n){return n.toLocaleString()}
function val(n){return n<0?'NA':n.toFixed(2)}
function esc(s){return s.replace(/[&<>"']/g,function(c){return{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]})}
function sourceInterval(n,begin,end){if(n.rv){const old=begin;begin=n.b-end;end=n.b-old}return esc(n.n)+'@'+fmt(begin)+'-'+fmt(end)+(n.rv?'-':'+')}
function tipText(e){return byId.get(e.l).n+' → '+byId.get(e.r).n+'\nvia '+byId.get(e.b).n+'\nphase '+val(e.pl)+' / '+val(e.prh)+' · support '+e.ss+'/'+e.sp+' · confidence '+e.cf.toFixed(2)}
function bindTip(mark,e){mark.onpointermove=function(event){tooltip.style.display='block';tooltip.style.left=Math.min(innerWidth-430,event.clientX+14)+'px';tooltip.style.top=Math.min(innerHeight-80,event.clientY+14)+'px';tooltip.textContent=tipText(e)};mark.onpointerleave=function(){tooltip.style.display='none'}}
function describe(e){const l=byId.get(e.l),r=byId.get(e.r),b=byId.get(e.b);return '<div><b>Status:</b> '+esc(e.st)+' &nbsp; <b>Component:</b> '+e.c+'</div><div class="coordinates"><b>Left:</b> '+sourceInterval(l,0,e.lp)+'<br><b>Right:</b> '+sourceInterval(r,e.rp,r.b)+'<br><b>Bridge:</b> '+sourceInterval(b,e.bl,e.br)+'</div><table class="evidence"><thead><tr><th>Boundary</th><th>Phase</th><th>Bubble bp</th><th>Minimizer</th><th>Identity</th><th>Coverage</th></tr></thead><tbody><tr><th>Left</th><td>'+val(e.pl)+'</td><td>'+fmt(e.plt)+'/'+fmt(e.plb)+'</td><td>'+val(e.ml)+'</td><td>'+val(e.il)+'</td><td>'+val(e.cl)+'</td></tr><tr><th>Right</th><td>'+val(e.prh)+'</td><td>'+fmt(e.prt)+'/'+fmt(e.prb)+'</td><td>'+val(e.mr)+'</td><td>'+val(e.ir)+'</td><td>'+val(e.cr)+'</td></tr></tbody></table><div class="support">Other hap spans gap: '+(e.hs?'yes':'no')+'<br>Sample support: '+e.ss+' / '+e.sp+'<br>Evidence: sample '+e.sscore.toFixed(2)+' · phase '+e.qscore.toFixed(2)+' · alignment '+e.ascore.toFixed(2)+'<br><b>Confidence: '+e.cf.toFixed(2)+'</b></div>'}
let selectedEdge=null;
function centerEdge(id){const marks=svg.querySelectorAll('[data-e="'+id+'"]');let left=Infinity,right=-Infinity,top=Infinity,bottom=-Infinity;marks.forEach(function(mark){const box=mark.getBBox();left=Math.min(left,box.x);right=Math.max(right,box.x+box.width);top=Math.min(top,box.y);bottom=Math.max(bottom,box.y+box.height)});if(left!==Infinity)view.scrollTo({left:Math.max(0,(left+right-view.clientWidth)/2),top:Math.max(0,(top+bottom-view.clientHeight)/2),behavior:'smooth'})}
function selectEdge(id){selectedEdge=selectedEdge===id?null:id;render();if(selectedEdge!==null)requestAnimationFrame(function(){centerEdge(id)})}
function anchorPosition(n,offset){return n.x+offset}
function anchorX(n,offset,scale,left,origin){return left+(anchorPosition(n,offset)-origin)*scale}
function render(){
  const cid=Number(comp.value),allNodes=N.filter(function(n){return n.c===cid}),allEdges=E.filter(function(e){return e.c===cid&&visible(e)}),focus=selectedEdge===null?null:E[selectedEdge];
  if(focus&&(focus.c!==cid||!allEdges.some(function(e){return e.i===focus.i}))){selectedEdge=null;return render()}
  const focused=focus!==null,nodeIds=focused?new Set([focus.l,focus.r,focus.b]):null;
  const nodes=focused?allNodes.filter(function(n){return nodeIds.has(n.i)}):allNodes;
  const edges=focused?[focus]:allEdges;
  const rows=[...new Set(nodes.map(function(n){return n.s+'\t'+n.h}))].sort(function(a,b){const x=a.split('\t'),y=b.split('\t');return x[0].localeCompare(y[0])||Number(x[1])-Number(y[1])});
  let origin=0,end=Math.max(1,...nodes.map(function(n){return n.x+n.b}));
  if(focused){const positions=[anchorPosition(byId.get(focus.l),focus.lp),anchorPosition(byId.get(focus.b),focus.bl),anchorPosition(byId.get(focus.b),focus.br),anchorPosition(byId.get(focus.r),focus.rp)],range=Math.max(1,Math.max(...positions)-Math.min(...positions)),padding=Math.max(500000,range*.16);origin=Math.max(0,Math.min(...positions)-padding);end=Math.max(...positions)+padding}
  const left=170,rowHeight=focused?86:58,top=58,span=Math.max(1,end-origin),factor=Number(zoom.value);
  const plotWidth=Math.max(1000,$('view').clientWidth-left-40)*factor,scale=plotWidth/span,width=left+plotWidth+35,height=Math.max(300,top+rows.length*rowHeight+45),pos=new Map();
  svg.replaceChildren();svg.setAttribute('viewBox','0 0 '+width+' '+height);svg.setAttribute('width',width);svg.setAttribute('height',height);zoomValue.textContent=factor+'x';
  rows.forEach(function(row,index){const parts=row.split('\t'),y=top+index*rowHeight;svg.appendChild(S('rect',{x:0,y:y-rowHeight/2,width:width,height:rowHeight,class:'lane '+(index%2?'alt':'')}));svg.appendChild(S('line',{x1:left,y1:y,x2:width-20,y2:y,class:'row'}));svg.appendChild(S('text',{x:14,y:y+4,class:'sample'},parts[0]+'  hap'+parts[1]));pos.set(row,y)});
  for(let tick=0;tick<=5;++tick){const bp=origin+span*tick/5,x=left+plotWidth*tick/5;svg.appendChild(S('line',{x1:x,y1:25,x2:x,y2:height-20,class:'row'}));svg.appendChild(S('text',{x:x,y:18,'text-anchor':'middle',class:'axis'},(bp/1000000).toFixed(span>=10000000?0:1)+' Mb'))}
  edges.forEach(function(e){const a=byId.get(e.l),m=byId.get(e.b),z=byId.get(e.r),active=focused?' active':'';if(!a||!m||!z)return;const points=[{x:anchorX(a,e.lp,scale,left,origin),y:pos.get(a.s+'\t'+a.h)},{x:anchorX(m,e.bl,scale,left,origin),y:pos.get(m.s+'\t'+m.h)},{x:anchorX(m,e.br,scale,left,origin),y:pos.get(m.s+'\t'+m.h)},{x:anchorX(z,e.rp,scale,left,origin),y:pos.get(z.s+'\t'+z.h)}];[[points[0],points[1]],[points[2],points[3]]].forEach(function(pair){const middle=(pair[0].y+pair[1].y)/2,path=S('path',{d:'M'+pair[0].x+','+pair[0].y+' C'+pair[0].x+','+middle+' '+pair[1].x+','+middle+' '+pair[1].x+','+pair[1].y,class:'edge '+(e.st==='keep'?'keep':'')+active,'data-e':e.i});path.onclick=function(){selectEdge(e.i)};bindTip(path,e);svg.appendChild(path)})});
  nodes.sort(function(a,b){return a.s.localeCompare(b.s)||a.h-b.h||a.x-b.x}).forEach(function(n){const y=pos.get(n.s+'\t'+n.h),nodeBegin=Math.max(n.x,origin),nodeEnd=Math.min(n.x+n.b,origin+span);if(nodeEnd<=nodeBegin)return;const x=left+(nodeBegin-origin)*scale,w=Math.max(1,(nodeEnd-nodeBegin)*scale),g=S('g',{class:'node h'+n.h+(focused?' active':''),'data-n':n.i});g.appendChild(S('rect',{x:x,y:y-(focused?16:11),width:w,height:focused?32:22,rx:3}));if(w>55)g.appendChild(S('text',{x:x+w/2,y:y+4,'text-anchor':'middle'},n.n.length>34?n.n.slice(0,31)+'...':n.n));g.appendChild(S('title',{},n.n+'\ncomponent: '+(n.x/1000000).toFixed(3)+'-'+((n.x+n.b)/1000000).toFixed(3)+' Mb\nlength: '+fmt(n.b)+' bp\norientation: '+(n.rv?'reverse':'forward')+'\n'+n.a+' -> '+n.z));g.onclick=function(event){event.stopPropagation();const edge=edges.find(function(e){return e.l===n.i||e.r===n.i||e.b===n.i});if(edge)selectEdge(edge.i)};svg.appendChild(g)});
  edges.forEach(function(e){const m=byId.get(e.b),y=pos.get(m.s+'\t'+m.h),x1=anchorX(m,e.bl,scale,left,origin),x2=anchorX(m,e.br,scale,left,origin),kind=(e.st==='keep'?'keep ':'')+(focused?'active':'');const bridgeMark=S('line',{x1:x1,y1:y,x2:x2,y2:y,class:'edge bridge-span '+kind,'data-e':e.i});bridgeMark.onclick=function(){selectEdge(e.i)};bindTip(bridgeMark,e);svg.appendChild(bridgeMark);[x1,x2].forEach(function(x){const dot=S('circle',{cx:x,cy:y,r:focused?6:4,class:'edge anchor '+kind,'data-e':e.i});dot.onclick=function(){selectEdge(e.i)};bindTip(dot,e);svg.appendChild(dot)})});
  list.replaceChildren();edges.slice().sort(function(a,b){return(a.st==='keep'?0:1)-(b.st==='keep'?0:1)||byId.get(a.l).x-byId.get(b.l).x}).forEach(function(e){const button=document.createElement('button'),route=document.createElement('span'),via=document.createElement('span'),metrics=document.createElement('span');button.className='candidate '+(e.st==='keep'?'keep ':'')+(focused?'active':'');button.dataset.e=e.i;route.className='route';route.textContent=byId.get(e.l).n+' → '+byId.get(e.r).n;via.className='via';via.textContent='via '+byId.get(e.b).n;metrics.className='metrics';metrics.textContent=e.st+' · confidence '+e.cf.toFixed(2)+' · support '+e.ss+'/'+e.sp;button.append(route,via,metrics);button.onclick=function(){selectEdge(e.i)};list.appendChild(button)});
  selection.classList.toggle('focused',focused);if(focused)detail.innerHTML=describe(focus);else detail.textContent='Click a connection or candidate.';const kept=edges.filter(function(e){return e.st==='keep'}).length;summary.textContent=focused?'Focused connection · '+nodes.length+' contigs':nodes.length+' contigs · '+kept+' kept · '+(edges.length-kept)+' rejected';view.scrollTo({left:0,top:0})
}
[...new Set(N.map(function(n){return n.c}))].sort(function(a,b){return a-b}).forEach(function(c){const o=document.createElement('option');o.value=c;o.textContent='component '+c;comp.appendChild(o)});comp.onchange=function(){selectedEdge=null;render()};filter.onchange=function(){selectedEdge=null;render()};query.oninput=function(){selectedEdge=null;render()};$('fit').onclick=function(){zoom.value=1;render()};$('clear').onclick=function(){selectedEdge=null;render()};let zoomFrame=0;zoom.oninput=function(){cancelAnimationFrame(zoomFrame);zoomFrame=requestAnimationFrame(render)};svg.onclick=function(event){if(!event.target.closest('[data-e],[data-n]')){selectedEdge=null;render()}};document.onkeydown=function(event){if(event.key==='Escape')$('clear').click()};render();
</script></body></html>)HTML");
    log_stream() << "  - HTML: " << file << "\n\n";
}

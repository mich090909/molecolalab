/**
 * MolecolaLab v0.4 — educational Cahn–Ingold–Prelog (CIP) support.
 * Supported: noncyclic structures; 3D tetrahedral centers C, Si, N+ with four
 * explicit single-bonded ligands, standard atomic numbers, elementary duplicate
 * node expansions for multiple bonds in the substituents. All unsupported
 * features are reported explicitly, not assigned an R/S descriptor.
 *
 * NOT a complete implementation of IUPAC CIP sequence rules (2013/2021).
 * Not supported: isotopes, cyclic ligands, pseudoasymmetric centers, E/Z-driven
 * ranking, descriptor-dependent stereogenic units, implicit H, or resonance.
 * 
 * This module does not import or distribute RDKit/Avogadro/Qt.
 */
import {element,findAtom} from './core.js';

const cipZ = a => element(a.symbol)?.Z ?? 0;
const cipPoint = a => [a.x,a.y,a.z];
const cipMinus = (a,b) => a.map((v,i)=>v-b[i]);
const cipCross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const cipDot=(a,b)=>a.reduce((t,v,i)=>t+v*b[i],0);
const cipNorm=a=>Math.hypot(...a);
const cipNeighbors=(graph,id)=>graph.get(id)||[];
function cipGraph(doc){
 const graph=new Map(doc.atoms.map(a=>[a.id,[]]));
 for(const b of doc.bonds){
  if(!graph.has(b.a)||!graph.has(b.b)||!Number.isInteger(b.order)||b.order<1||b.order>3)return null;
  graph.get(b.a).push({id:b.b,order:b.order});
  graph.get(b.b).push({id:b.a,order:b.order});
 }
 return graph;
}
function cipHasCycle(graph,from){
 const seen=new Set(),stack=[[from,null]];
 while(stack.length){const [id,parent]=stack.pop();if(seen.has(id))return true;seen.add(id);
  for(const other of cipNeighbors(graph,id)){if(other.id===parent)continue;
   if(seen.has(other.id))return true;stack.push([other.id,id]);
  }
 }
 return false;
}
/** Atomic-number levels from one ligand, with duplicate terminal nodes for C=C/C=O/C≡N.
 * Valid on acyclic explicit molecular graphs only.
 */
function cipLayers(doc,graph,centerId,ligandId){
 const atoms=new Map(doc.atoms.map(a=>[a.id,a]));
 let frontier=[{id:ligandId,parent:centerId,ghost:false}];
 const levels=[];
 for(let depth=0;depth<=doc.atoms.length+1&&frontier.length;depth++){
  const labels=frontier.map(p=>cipZ(atoms.get(p.id))).sort((a,b)=>b-a);
  levels.push(labels);
  const next=[];
  for(const node of frontier){
   if(node.ghost)continue;
   const edges=cipNeighbors(graph,node.id);
   // The back-duplicate of the parent on a multiple bond is terminal in CIP's digraph.
   const back=edges.find(e=>e.id===node.parent);
   if(back && back.order>1)for(let j=1;j<back.order;j++)next.push({id:node.parent,parent:node.id,ghost:true});
   for(const e of edges){
    if(e.id===node.parent)continue;
    next.push({id:e.id,parent:node.id,ghost:false});
    for(let j=1;j<e.order;j++)next.push({id:e.id,parent:node.id,ghost:true});
   }
  }
  frontier=next;
  if(frontier.length>100000)throw Error('Confronto CIP troppo esteso');
 }
 return levels;
}
function cipCompareLayers(one,two){
 for(let depth=0;depth<Math.max(one.length,two.length);depth++){
  const a=one[depth]||[],b=two[depth]||[];
  for(let i=0;i<Math.max(a.length,b.length);i++){
   const diff=(a[i]||0)-(b[i]||0);
   if(diff!==0)return {order:Math.sign(diff),depth};
  }
 }
 return {order:0,depth:null};
}
/** Identify chemically identical ligand trees in the supported acyclic subset.
 * A CIP-level tie that is not a trivial structural duplication is unresolved.
 */
function cipBranchSignature(doc,graph,node,parent){
 const a=findAtom(doc,node),parts=[];
 for(const b of cipNeighbors(graph,node))if(b.id!==parent)
  parts.push(`${b.order}:${cipBranchSignature(doc,graph,b.id,node)}`);
 return `${a.symbol}:${Number(a.charge||0)}(${parts.sort().join('|')})`;
}
/** Main educational API. `assigned` is never emitted for an unsupported graph. */
export function analyzeCIP(doc,id){
 const central=findAtom(doc,id);
 if(!central)return {status:'unavailable',reason:'Atomo non trovato.'};
 const graph=cipGraph(doc);
 if(!graph)return {status:'unsupported',reason:'Connettività non valida.'};
 const adjacent=cipNeighbors(graph,id);
 if(adjacent.length!==4||adjacent.some(e=>e.order!==1))
  return {status:'not_tetrahedral',reason:'Occorrono quattro sostituenti espliciti, uniti mediante quattro legami singoli.'};
 if(!['C','Si'].includes(central.symbol)&&!(central.symbol==='N'&&Number(central.charge||0)===1))
  return {status:'unsupported',reason:'Assegnazione v0.4 limitata a C, Si e N⁺ tetraedrici con quattro sostituenti.'};
 if(['C','Si'].includes(central.symbol)&&Number(central.charge||0)!==0)
  return {status:'unsupported',reason:'Centro C/Si formalmente carico: non supportato in questa versione.'};
 if(doc.atoms.some(a=>a.isotope!=null||a.massNumber!=null||a.mass!=null))
  return {status:'unsupported',reason:'Priorità isotopiche non implementate.'};
 if(doc.atoms.length>250)
  return {status:'unsupported',reason:'Analisi CIP limitata a 250 atomi.'};
 if(cipHasCycle(graph,id))
  return {status:'unsupported',reason:'La molecola contiene un anello: il confronto CIP dei percorsi ciclici non è ancora implementato.'};
 const entries=adjacent.map(({id:neighbor})=>{
  const layers=cipLayers(doc,graph,id,neighbor);
  const atom=findAtom(doc,neighbor);
  return {id:neighbor,symbol:atom.symbol,layers};
 });
 entries.sort((a,b)=>cipCompareLayers(b.layers,a.layers).order || a.id-b.id);
 const ties=[];
 for(let i=0;i<3;i++)if(cipCompareLayers(entries[i].layers,entries[i+1].layers).order===0)ties.push([entries[i].id,entries[i+1].id]);
 if(ties.length){
  const equivalent=ties.every(([a,b])=>cipBranchSignature(doc,graph,a,id)===cipBranchSignature(doc,graph,b,id));
  // Do not declare true equivalence in the presence of stereogenic units elsewhere:
  // enantiomorphic ligands can require descriptor-dependent CIP sequence rules.
  const otherPotential=ties.some(([left,right])=>[left,right].some(first=>{
   const stack=[[first,id]];
   while(stack.length){const [node,parent]=stack.pop();const current=findAtom(doc,node);
    if(cipNeighbors(graph,node).length===4 &&
     (['C','Si'].includes(current.symbol)||(current.symbol==='N'&&current.charge===1)) &&
     new Set(cipNeighbors(graph,node).map(e=>findAtom(doc,e.id).symbol)).size>1)return true;
    for(const e of cipNeighbors(graph,node))if(e.id!==parent)stack.push([e.id,node]);
   }
   return false;
  }));
  return equivalent&&!otherPotential
   ? {status:'not_stereogenic',reason:'Sostituenti equivalenti: nessun centro stereogenico R/S (nel modello supportato).'}
   : {status:'unresolved',reason:'Priorità indistinguibili con le regole attualmente disponibili (possibili condizioni stereochimiche complesse): nessuna assegnazione R/S.',ties};
 }
 const points=entries.map(e=>cipPoint(findAtom(doc,e.id)));
 if(points.some(p=>p.some(v=>!Number.isFinite(v))))return {status:'unsupported',reason:'Coordinate non valide.'};
 const a=cipMinus(points[0],points[3]),b=cipMinus(points[1],points[3]),c=cipMinus(points[2],points[3]);
 const denom=cipNorm(a)*cipNorm(b)*cipNorm(c);
 const orientation=cipDot(a,cipCross(b,c));
 if(!Number.isFinite(orientation)||denom<1e-10||Math.abs(orientation/denom)<0.06)
  return {status:'degenerate',reason:'Geometria tridimensionale planare o quasi degenerata: impossibile assegnare R/S in modo attendibile.',priorities:entries.map((v,i)=>({rank:i+1,atomId:v.id,symbol:v.symbol}))};
 return {status:'assigned',descriptor:orientation<0?'R':'S',centerId:id,symbol:central.symbol,
  priorities:entries.map((v,i)=>({rank:i+1,atomId:v.id,symbol:v.symbol})),
  geometryScore:Math.abs(orientation/denom),
  note:'Assegnazione CIP didattica basata sul grafo aciclico e sulle coordinate 3D; tutti e quattro i sostituenti devono essere espliciti. Non sono gestiti isotopi, anelli e regole CIP avanzate.'};
}
/** Only atoms which meet the local connectivity requirement are listed, including unresolved cases. */
export function summarizeCIP(doc){
 const g=cipGraph(doc);
 if(!g)return [];
 return doc.atoms.filter(a=>cipNeighbors(g,a.id).length===4&&cipNeighbors(g,a.id).every(e=>e.order===1)&&(['C','Si'].includes(a.symbol)||(a.symbol==='N'&&a.charge===1)))
  .map(a=>({id:a.id,symbol:a.symbol,...analyzeCIP(doc,a.id)}))
  .filter(a=>a.status!=='not_stereogenic');
}
/** Spatial inversion. Mirrors all atomic coordinates, not just one stereocenter. */
export function mirrorMolecule(doc){
 if(!doc.atoms.length)throw Error('Non ci sono atomi da riflettere.');
 const x0=doc.atoms.reduce((s,a)=>s+a.x,0)/doc.atoms.length;
 for(const a of doc.atoms)a.x=2*x0-a.x;
 doc.name='Speculare di '+(doc.name||'molecola').replace(/^Speculare di /,'');
}

/**
 * MolecolaLab v0.5 — E/Z for explicitly modelled acyclic C=C double bonds.
 * Uses the same atomic-number first CIP ligand comparison as R/S. This is NOT
 * the full 2013/2021 CIP rule set: isotope, ring, cumulene, implicit H,
 * stereodescriptor-dependent tie-breaks and exotic bonding are unresolved.
 * Geometry is determined from 3D coordinates, NOT from screen projection.
 */
const ezSubtract=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]];
const ezUnit=v=>{const n=cipNorm(v);return n>1e-9?v.map(x=>x/n):null;};
const ezPerpendicular=(point,origin,axis)=>{
 const p=ezSubtract(cipPoint(point),cipPoint(origin));
 const axial=cipDot(p,axis);
 return ezUnit(p.map((x,i)=>x-axial*axis[i]));
};
const ezCharge=a=>Number(a.charge||0);

/** Assign E or Z, or explicitly report why assignment is not justified. */
export function analyzeEZ(doc,bondId){
 const bond=doc.bonds.find(b=>b.id===bondId);
 if(!bond)return {status:'unavailable',reason:'Legame non trovato.'};
 if(bond.order!==2)return {status:'not_double',reason:'Seleziona un doppio legame C=C.'};
 const first=findAtom(doc,bond.a),second=findAtom(doc,bond.b);
 if(!first||!second)return {status:'unsupported',reason:'Estremi del doppio legame mancanti.'};
 if(first.symbol!=='C'||second.symbol!=='C')return {status:'unsupported',reason:'La versione 0.5 assegna E/Z ai doppi legami carbonio–carbonio (C=C) soltanto.'};
 if(doc.atoms.length>250)return {status:'unsupported',reason:'Analisi E/Z limitata a 250 atomi.'};
 if(doc.atoms.some(a=>a.isotope!=null||a.massNumber!=null||a.mass!=null))
  return {status:'unsupported',reason:'Priorità di isotopi non supportate.'};
 const graph=cipGraph(doc);
 if(!graph)return {status:'unsupported',reason:'Connettività non valida.'};
 if(cipHasCycle(graph,bond.a))return {status:'unsupported',reason:'Molecola con anelli: confronto CIP ciclico non supportato.'};
 if(doc.atoms.some(a=>ezCharge(a)!==0))
  return {status:'unsupported',reason:'La procedura E/Z didattica non gestisce cariche formali.'};
 const ends=[first,second];
 const priorityEnds=[];
 for(const end of ends){
  const neighbors=cipNeighbors(graph,end.id),side=neighbors.filter(e=>e.id!==(end.id===first.id?second.id:first.id));
  if(neighbors.length!==3||side.length!==2||neighbors.reduce((s,e)=>s+e.order,0)!==4||side.some(e=>e.order!==1))
   return {status:'unsupported',reason:'Servono due sostituenti espliciti con legami singoli su ciascun carbonio sp² (anche gli H vanno disegnati).'};
  const ligands=side.map(e=>{
   const at=findAtom(doc,e.id);
   return {id:e.id,symbol:at.symbol,layers:cipLayers(doc,graph,end.id,e.id)};
  });
  ligands.sort((a,b)=>cipCompareLayers(b.layers,a.layers).order||a.id-b.id);
  const comparison=cipCompareLayers(ligands[0].layers,ligands[1].layers).order;
  if(comparison===0){
   const identical=cipBranchSignature(doc,graph,ligands[0].id,end.id)===cipBranchSignature(doc,graph,ligands[1].id,end.id);
   return identical
    ? {status:'not_stereogenic',reason:`I due sostituenti del carbonio ${end.id} sono uguali: E/Z non definito.`}
    : {status:'unresolved',reason:`Non riesco a distinguere con certezza le priorità CIP sul carbonio ${end.id}. Non assegno E/Z.`};
  }
  priorityEnds.push({endId:end.id,ranked:ligands.map((a,i)=>({rank:i+1,atomId:a.id,symbol:a.symbol}))});
 }
 const axis=ezUnit(ezSubtract(cipPoint(second),cipPoint(first)));
 if(!axis)return {status:'degenerate',reason:'I due carboni del doppio legame coincidono nelle coordinate.'};
 const oriented=priorityEnds.map((entry,i)=>{
  const base=ends[i];
  return entry.ranked.map(s=>ezPerpendicular(findAtom(doc,s.atomId),base,axis));
 });
 if(oriented.some(pair=>pair.some(v=>!v)))
  return {status:'degenerate',reason:'Uno o più sostituenti sono allineati al doppio legame: geometria non affidabile.',ends:priorityEnds};
 if(oriented.some(pair=>cipDot(pair[0],pair[1])>-0.45))
  return {status:'degenerate',reason:'Geometria locale non compatibile con i due sostituenti sp² su lati opposti: correggi le coordinate.',ends:priorityEnds};
 const dot=cipDot(oriented[0][0],oriented[1][0]);
 if(!Number.isFinite(dot)||Math.abs(dot)<0.85)
  return {status:'degenerate',reason:'Doppio legame molto deformato o sostituenti fuori piano: nessuna assegnazione E/Z attendibile.',ends:priorityEnds};
 return {status:'assigned',descriptor:dot>0?'Z':'E',bondId:bond.id,firstId:first.id,secondId:second.id,
  ends:priorityEnds,geometryScore:Math.abs(dot),
  note:'I sostituenti di priorità CIP 1 si trovano dallo stesso lato (Z) oppure su lati opposti (E) del legame C=C. Calcolo limitato a grafi aciclici con tutti gli H espliciti e geometria sp² attendibile.'};
}

/** All C=C bonds for the E/Z panel, including ties and unresolved geometry. */
export function summarizeEZ(doc){
 return doc.bonds.filter(b=>b.order===2&&findAtom(doc,b.a)?.symbol==='C'&&findAtom(doc,b.b)?.symbol==='C')
  .map(b=>({bondId:b.id,firstId:b.a,secondId:b.b,...analyzeEZ(doc,b.id)}));
}

/**
 * Construct the OPPOSITE geometric isomer by rotating the COMPLETE component
 * on one side of a formally disconnected alkene bond by 180° around the C=C axis.
 * This is a rebuild of a different isomer, NOT a physically allowed internal
 * torsional rotation around the double bond. Atomic positions and connections
 * on each side are preserved rigidly.
 */
export function invertEZ(doc,bondId){
 const before=analyzeEZ(doc,bondId);
 if(before.status!=='assigned')throw Error('E/Z non assegnabile: '+before.reason);
 const bond=doc.bonds.find(b=>b.id===bondId),base=findAtom(doc,bond.a),pivot=findAtom(doc,bond.b);
 const axis=ezUnit(ezSubtract(cipPoint(pivot),cipPoint(base)));
 if(!axis)throw Error('Asse C=C degenerato.');
 const graph=cipGraph(doc),seen=new Set([pivot.id]),stack=[pivot.id];
 while(stack.length){const now=stack.pop();for(const entry of cipNeighbors(graph,now)){
   // The C=C is disconnected for this traversal.
   if((now===pivot.id&&entry.id===base.id)||(now===base.id&&entry.id===pivot.id))continue;
   if(!seen.has(entry.id)){seen.add(entry.id);stack.push(entry.id);}
 }}
 if(seen.has(base.id))throw Error('Doppio legame all’interno di un anello: impossibile generare l’isomero opposto.');
 const coords=new Map();
 for(const id of seen){const atom=findAtom(doc,id),diff=ezSubtract(cipPoint(atom),cipPoint(pivot));
  const proj=cipDot(diff,axis),rot=diff.map((v,i)=>2*proj*axis[i]-v);
  coords.set(id,rot.map((v,i)=>v+cipPoint(pivot)[i]));
 }
 const original=new Map([...seen].map(id=>{const a=findAtom(doc,id);return [id,[a.x,a.y,a.z]];}));
 for(const [id,xyz] of coords){const a=findAtom(doc,id);[a.x,a.y,a.z]=xyz;}
 const after=analyzeEZ(doc,bondId);
 if(after.status!=='assigned'||after.descriptor===before.descriptor){
  for(const [id,xyz] of original){const a=findAtom(doc,id);[a.x,a.y,a.z]=xyz;}
  throw Error('Impossibile costruire con affidabilità l’altro isomero geometrico.');
 }
 if(/^\([EZ]\)-/.test(doc.name||''))doc.name=`(${after.descriptor})-${doc.name.slice(4)}`;
 else doc.name=`Isomero ${after.descriptor} di ${doc.name||'molecola'}`.slice(0,80);
 return {before:before.descriptor,after:after.descriptor,bondId,rotatedAtoms:seen.size-1};
}

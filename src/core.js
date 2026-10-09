/** MolecolaLab Web: data model and numerical operations. Independent educational prototype. */
export const SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split(' ');
if(SYMBOLS.length !== 118) throw new Error('Errore tavola periodica');
const data = {
 H:[0.31,1,'#f9fafb'], He:[0.28,0,'#f5bde6'], Li:[1.28,1,'#c29dff'], Be:[0.96,2,'#8cddb9'],
 B:[0.84,3,'#f8b879'], C:[0.76,4,'#303d52'], N:[0.71,3,'#4287ed'], O:[0.66,2,'#eb5261'],
 F:[0.57,1,'#5dc99d'], Ne:[0.58,0,'#8ed9e8'], Na:[1.66,1,'#b88dfa'], Mg:[1.41,2,'#7fce9d'],
 Al:[1.21,3,'#ced3de'], Si:[1.11,4,'#d19d71'], P:[1.07,3,'#ec9b3e'], S:[1.05,2,'#edc849'],
 Cl:[1.02,1,'#51bd76'], Ar:[1.06,0,'#7ad0ef'], K:[2.03,1,'#a96fde'], Ca:[1.76,2,'#a1d4b2'],
 Fe:[1.24,null,'#b99581'], Cu:[1.32,null,'#bd7c55'], Zn:[1.22,null,'#a7b7cd'],
 Br:[1.20,1,'#a65b43'], I:[1.39,1,'#9364b9'], Xe:[1.40,0,'#7ea7cf'],
};
export const ELEMENTS=SYMBOLS.map((symbol,i)=>({symbol,Z:i+1,radius:data[symbol]?.[0]??1.24,valence:data[symbol]?.[1]??null,color:data[symbol]?.[2]??'#b7c1d1'}));
export function element(symbol){return ELEMENTS.find(e=>e.symbol===symbol)||null;}
export function atomicElement(Z){return ELEMENTS[Number(Z)-1]||null;}
export function blank(){return {atoms:[],bonds:[],name:'Nuova molecola'};}
export function clone(doc){return JSON.parse(JSON.stringify(doc));}
export function atomId(doc){return Math.max(0,...doc.atoms.map(a=>a.id))+1;}
export function bondId(doc){return Math.max(0,...doc.bonds.map(b=>b.id))+1;}
export function findAtom(doc,id){return doc.atoms.find(a=>a.id===id);}
export function bondBetween(doc,a,b){return doc.bonds.find(k=>(k.a===a&&k.b===b)||(k.a===b&&k.b===a));}
export function addAtom(doc,symbol,x,y,z,connectedTo=null,bondOrder=1){
 if(!element(symbol))throw new Error('Elemento sconosciuto: '+symbol);
 if(![x,y,z].every(Number.isFinite))throw new Error('Coordinate non valide');
 if(connectedTo!=null&&(!findAtom(doc,connectedTo)||![1,2,3].includes(bondOrder)))throw new Error('Collegamento o ordine di legame non valido');
 const id=atomId(doc);doc.atoms.push({id,symbol,x,y,z,charge:0});
 if(connectedTo!=null)doc.bonds.push({id:bondId(doc),a:connectedTo,b:id,order:bondOrder});
 return id;
}
export function setBond(doc,a,b,order=1){
 if(a===b||!findAtom(doc,a)||!findAtom(doc,b))throw new Error('Seleziona due atomi distinti');
 if(![1,2,3].includes(order))throw new Error('Ordine di legame non supportato');
 const existing=bondBetween(doc,a,b);if(existing) existing.order=order;else doc.bonds.push({id:bondId(doc),a,b,order});
}
export function eraseAtom(doc,id){doc.atoms=doc.atoms.filter(a=>a.id!==id);doc.bonds=doc.bonds.filter(b=>b.a!==id&&b.b!==id);}
export function eraseBond(doc,id){doc.bonds=doc.bonds.filter(b=>b.id!==id);}
/** Geometry edits operate on coordinates and never imply an energy minimization. */
export function setAtomPosition(doc,id,x,y,z){
 const atom=findAtom(doc,id);
 if(!atom)throw Error('Atomo inesistente');
 if(![x,y,z].every(n=>Number.isFinite(n)&&Math.abs(n)<=1000))throw Error('Coordinate ammesse da −1000 a +1000 Å');
 Object.assign(atom,{x,y,z});
}
export function setAtomElement(doc,id,symbol){
 const a=findAtom(doc,id);if(!a||!element(symbol))throw Error('Atomo o elemento non valido');
 a.symbol=symbol;
}
export function setBondLength(doc,id,length){
 const b=doc.bonds.find(b=>b.id===id);
 if(!b)throw Error('Legame non trovato');
 if(!Number.isFinite(length)||length<.2||length>6)throw Error('Lunghezza richiesta tra 0,20 e 6 Å');
 const a=findAtom(doc,b.a),c=findAtom(doc,b.b),d=dist(a,c);
 if(d<1e-8)throw Error('Atomi coincidenti: sposta un atomo prima di modificare il legame');
 setAtomPosition(doc,c.id,a.x+(c.x-a.x)*length/d,a.y+(c.y-a.y)*length/d,a.z+(c.z-a.z)*length/d);
}
export const netCharge=doc=>doc.atoms.reduce((s,a)=>s+(Number(a.charge)||0),0);
/** Only simple main-group Lewis structures; experimental geometry comes from coordinate measurements. */
export function estimateVSEPR(doc,id){
 const atom=findAtom(doc,id);
 if(!atom)return null;
 const neighbors=doc.bonds.filter(b=>b.a===id||b.b===id);
 const n=neighbors.length;
 if(n<2)return {supported:false,reason:'Servono almeno due atomi legati al centro per identificare una geometria.'};
 const ve={B:3,C:4,N:5,O:6,F:7,Si:4,P:5,S:6,Cl:7,Br:7,I:7,Xe:8};
 if(!Object.hasOwn(ve,atom.symbol))return {supported:false,reason:'Elemento non incluso nel modello VSEPR semplificato.'};
 const charge=Number(atom.charge||0);
 const bondedElectrons=neighbors.reduce((s,b)=>s+b.order,0);
 const remaining=ve[atom.symbol]-charge-bondedElectrons;
 if(remaining<0||remaining%2!==0)return {supported:false,reason:'Elettroni non appaiati o conteggio Lewis non compatibile: verificare legami e cariche.'};
 const lonePairs=remaining/2,domains=n+lonePairs;
 const electronic={2:'lineare',3:'trigonale planare',4:'tetraedrica',5:'bipiramidale trigonale',6:'ottaedrica'};
 const shapes={
  '2:0':'lineare','2:1':'angolare','2:2':'angolare','2:3':'lineare',
  '3:0':'trigonale planare','3:1':'piramidale trigonale','3:2':'a T', '3:3':'a T',
  '4:0':'tetraedrica','4:1':'a cavalletto (altalena)','4:2':'quadrata planare',
  '5:0':'bipiramidale trigonale','5:1':'piramidale quadrata',
  '6:0':'ottaedrica'
 };
 const shape=shapes[`${n}:${lonePairs}`];
 if(!electronic[domains]||!shape)return {supported:false,reason:'Geometria fuori dal campo AXₙEₘ supportato.'};
 return {supported:true,central:atom.symbol,bondedAtoms:n,lonePairs,electronDomains:domains,formula:`AX${n}${lonePairs?'E'+(lonePairs===1?'':lonePairs):''}`,
  electronicGeometry:electronic[domains],molecularGeometry:shape,
  note:'Stima Lewis/VSEPR per specie semplici a guscio chiuso; non sostituisce un calcolo elettronico né verifica risonanza o stabilità.'};
}

export const dist=(a,b)=>Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z);
const sub=(a,b)=>[a.x-b.x,a.y-b.y,a.z-b.z];
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const mag=a=>Math.hypot(...a);
const unit=a=>{const l=mag(a)||1;return a.map(v=>v/l)};
export function angle(a,b,c){const u=sub(a,b),v=sub(c,b);if(mag(u)<1e-8||mag(v)<1e-8)return NaN;return Math.acos(Math.max(-1,Math.min(1,dot(u,v)/mag(u)/mag(v))))*180/Math.PI;}
export function dihedral(p1,p2,p3,p4){
 const b1=sub(p2,p1),b2=sub(p3,p2),b3=sub(p4,p3);
 const n1=cross(b1,b2),n2=cross(b2,b3);
 if(mag(n1)<1e-10||mag(n2)<1e-10)return NaN;
 const m1=cross(unit(n1),unit(b2));return Math.atan2(dot(m1,unit(n2)),dot(unit(n1),unit(n2)))*180/Math.PI;
}
export function formula(doc){const counts={};for(const a of doc.atoms)counts[a.symbol]=(counts[a.symbol]||0)+1;
 const keys=Object.keys(counts).sort((a,b)=>a==='C'?-1:b==='C'?1:a==='H'?-1:b==='H'?1:a.localeCompare(b));
 return keys.map(s=>s+(counts[s]>1?counts[s]:'' )).join('')||'—';}
export function centroid(doc){if(!doc.atoms.length)return {x:0,y:0,z:0};return doc.atoms.reduce((p,a)=>({x:p.x+a.x/doc.atoms.length,y:p.y+a.y/doc.atoms.length,z:p.z+a.z/doc.atoms.length}),{x:0,y:0,z:0});}
export function center(doc){const c=centroid(doc);for(const a of doc.atoms){a.x-=c.x;a.y-=c.y;a.z-=c.z;}}
export function suggestedPosition(doc,anchorId){
 const a=findAtom(doc,anchorId);if(!a)return {x:0,y:0,z:0};
 const neighbors=doc.bonds.filter(b=>b.a===a.id||b.b===a.id).map(b=>findAtom(doc,b.a===a.id?b.b:b.a));
 let v=[1,0,0];if(neighbors.length){const sums=neighbors.reduce((s,n)=>{const d=unit(sub(n,a));return s.map((k,i)=>k-d[i]);},[0,0,0]);v=mag(sums)>.15?unit(sums):unit([0.3,0.7,0.85]);}
 const seed=neighbors.length*1.17;v=unit([v[0]+0.19*Math.sin(seed),v[1]+0.23*Math.cos(seed),v[2]+0.25*Math.sin(seed+1)]);
 return {x:a.x+v[0]*1.46,y:a.y+v[1]*1.46,z:a.z+v[2]*1.46};
}
export function hydrogenate(doc){
 // Simplified valence completion: neutral main-group single-center bonding only.
 let total=0;
 for(const a of [...doc.atoms]){
  if(a.symbol==='H'||Number(a.charge||0)!==0)continue;
  const val=element(a.symbol)?.valence;if(val==null||val===0)continue;
  const neighbors=doc.bonds.filter(b=>b.a===a.id||b.b===a.id).map(b=>({b,other:findAtom(doc,b.a===a.id?b.b:b.a)}));
  const missing=Math.max(0,Math.min(4,Math.round(val-neighbors.reduce((s,n)=>s+n.b.order,0))));
  for(let i=0;i<missing;i++){
   const attached=doc.bonds.filter(b=>b.a===a.id||b.b===a.id).map(b=>findAtom(doc,b.a===a.id?b.b:b.a));
   let best=[1,0,0],score=-Infinity;
   // Deterministic spherical search for maximum angular separation.
   for(let k=0;k<110;k++){
    const z=1-2*(k+.5)/110,phi=k*2.399963229728653;
    const candidate=[Math.sqrt(1-z*z)*Math.cos(phi),Math.sqrt(1-z*z)*Math.sin(phi),z];
    let s=0;for(const other of attached){const u=unit(sub(other,a));s+=Math.pow(1-dot(candidate,u),2);}
    // Spread molecules deterministically if the atom has no neighbors.
    if(!attached.length)s=candidate[0]*.0001;
    if(s>score){score=s;best=candidate;}
   }
   const len=(element(a.symbol).radius+element('H').radius)*1.12;
   addAtom(doc,'H',a.x+len*best[0],a.y+len*best[1],a.z+len*best[2],a.id,1);total++;
  }
 }
 return total;
}
export function autoBond(doc){
 // XYZ cannot encode bonds; simple distance-based guess, not bond-order perception.
 let made=0;
 for(let i=0;i<doc.atoms.length;i++)for(let j=i+1;j<doc.atoms.length;j++){
  const a=doc.atoms[i],b=doc.atoms[j];const d=dist(a,b),max=(element(a.symbol).radius+element(b.symbol).radius)*1.25;
  if(d>.25&&d<=max&&!bondBetween(doc,a.id,b.id)){
   const n1=doc.bonds.filter(k=>k.a===a.id||k.b===a.id).length;
   const n2=doc.bonds.filter(k=>k.a===b.id||k.b===b.id).length;
   const v1=element(a.symbol).valence??8,v2=element(b.symbol).valence??8;
   if(n1<v1&&n2<v2){setBond(doc,a.id,b.id,1);made++;}
  }
 }
 return made;
}
export function regularize(doc,iterations=160){
 // Graphical spacing only; NOT a molecular mechanics force field.
 if(doc.atoms.length<2)return;
 const origin=centroid(doc), bonded=new Set(doc.bonds.map(b=>[Math.min(b.a,b.b),Math.max(b.a,b.b)].join(':')));
 for(let t=0;t<iterations;t++){
  const forces=new Map(doc.atoms.map(a=>[a.id,[0,0,0]]));
  for(const b of doc.bonds){const a=findAtom(doc,b.a),c=findAtom(doc,b.b);const delta=sub(c,a),d=mag(delta)||1e-5;
   const target=(element(a.symbol).radius+element(c.symbol).radius)*(b.order===3?.87:b.order===2?.91:1.05);
   const e=(d-target)*.035;const f=delta.map(v=>v/d*e);for(let k=0;k<3;k++){forces.get(a.id)[k]+=f[k];forces.get(c.id)[k]-=f[k];}
  }
  for(let i=0;i<doc.atoms.length;i++)for(let j=i+1;j<doc.atoms.length;j++){
   const a=doc.atoms[i],b=doc.atoms[j],key=[Math.min(a.id,b.id),Math.max(a.id,b.id)].join(':');if(bonded.has(key))continue;
   let delta=sub(b,a),d=mag(delta);
   if(d<.01){delta=[.1,(j%3+1)*.12,.08];d=mag(delta);}
   const min=(element(a.symbol).radius+element(b.symbol).radius)*1.35;
   if(d<min){const f=delta.map(v=>v/d*(min-d)*.014);for(let k=0;k<3;k++){forces.get(a.id)[k]-=f[k];forces.get(b.id)[k]+=f[k];}}
  }
  for(const a of doc.atoms){const f=forces.get(a.id);a.x+=f[0];a.y+=f[1];a.z+=f[2];}
 }
 const after=centroid(doc);for(const a of doc.atoms){a.x+=origin.x-after.x;a.y+=origin.y-after.y;a.z+=origin.z-after.z;}
}
const xyznum=(n)=>Number(n).toFixed(5);
export function exportXYZ(doc){return `${doc.atoms.length}\n${doc.name||'MolecolaLab'} — formato XYZ; la connettivita non e inclusa\n`+doc.atoms.map(a=>`${a.symbol} ${xyznum(a.x)} ${xyznum(a.y)} ${xyznum(a.z)}`).join('\n')+'\n';}
export function importXYZ(content){
 const lines=content.replace(/\r/g,'').split('\n');const n=Number(lines[0]?.trim());if(!Number.isInteger(n)||n<0||n>6000)throw Error('XYZ: conteggio atomi non valido');
 if(lines.length<n+2)throw Error('XYZ: file incompleto');const doc=blank();doc.name=(lines[1]||'Molecola importata').slice(0,80);
 for(let i=0;i<n;i++){
  const p=lines[i+2].trim().split(/\s+/);const sym=/^\d+$/.test(p[0])?atomicElement(+p[0])?.symbol:p[0];if(p.length<4||!element(sym))throw Error('XYZ: elemento non valido alla riga '+(i+3));
  const coords=p.slice(1,4).map(Number);if(!coords.every(Number.isFinite))throw Error('XYZ: coordinate non valide alla riga '+(i+3));
  addAtom(doc,sym,...coords);
 }
 autoBond(doc);return doc;
}
const pad=(n,w)=>String(n).padStart(w,' ');
const molcoord=n=>Number(n).toFixed(4).padStart(10,' ');
export function exportMOL(doc){
 if(doc.atoms.length>999||doc.bonds.length>999)throw Error('MOL V2000 ammette massimo 999 atomi e legami');
 const lines=[(doc.name||'MolecolaLab').slice(0,80),'  MolecolaLab Web v0.2','',`${pad(doc.atoms.length,3)}${pad(doc.bonds.length,3)}  0  0  0  0            999 V2000`];
 for(const a of doc.atoms)lines.push(`${molcoord(a.x)}${molcoord(a.y)}${molcoord(a.z)} ${a.symbol.padEnd(3,' ')} 0  0  0  0  0  0  0  0  0  0  0  0`);
 for(const b of doc.bonds)lines.push(`${pad(doc.atoms.findIndex(a=>a.id===b.a)+1,3)}${pad(doc.atoms.findIndex(a=>a.id===b.b)+1,3)}${pad(b.order,3)}  0  0  0  0`);
 const charged=doc.atoms.map((a,i)=>({i:i+1,charge:a.charge||0})).filter(c=>c.charge);
 for(let i=0;i<charged.length;i+=8){const slice=charged.slice(i,i+8);lines.push(`M  CHG${pad(slice.length,3)}${slice.map(c=>pad(c.i,4)+pad(c.charge,4)).join('')}`);}
 lines.push('M  END');return lines.join('\n')+'\n';
}
export function importMOL(content){
 const lines=content.replace(/\r/g,'').split('\n');if(lines.length<5)throw Error('MOL incompleto');
 const counts=lines[3];if(!/V2000/.test(counts))throw Error('Supportato solo MOL V2000');
 const n=Number(counts.slice(0,3)),m=Number(counts.slice(3,6));if(!Number.isInteger(n)||!Number.isInteger(m)||n<0||m<0||n>999||m>999)throw Error('Conteggio MOL non valido');
 if(lines.length<4+n+m)throw Error('MOL: righe insufficienti');
 const doc=blank();doc.name=(lines[0]||'Molecola MOL').trim().slice(0,80);
 for(let i=0;i<n;i++){
  const line=lines[4+i],sym=line.slice(31,34).trim(),coords=[line.slice(0,10),line.slice(10,20),line.slice(20,30)].map(v=>Number(v.trim()));
  if(!element(sym)||!coords.every(Number.isFinite))throw Error('MOL: atomo '+(i+1)+' non valido');addAtom(doc,sym,...coords);
 }
 for(let j=0;j<m;j++){
  const ln=lines[4+n+j],a=Number(ln.slice(0,3)),b=Number(ln.slice(3,6)),order=Number(ln.slice(6,9));
  if(a<1||a>n||b<1||b>n||!Number.isInteger(a)||!Number.isInteger(b)||![1,2,3].includes(order))throw Error('MOL: legame '+(j+1)+' non valido');setBond(doc,a,b,order);
 }
 for(const line of lines.slice(4+n+m))if(line.startsWith('M  CHG')){
  const pairs=Number(line.slice(6,9));for(let i=0;i<pairs;i++){const idx=Number(line.slice(9+i*8,13+i*8)),charge=Number(line.slice(13+i*8,17+i*8));if(findAtom(doc,idx))findAtom(doc,idx).charge=charge;}
 }
 return doc;
}
/** Chemical JSON (CJSON) interchange: original implementation following OpenChemistry's public schema.
 * Supported: atomic numbers, Cartesian coordinates, formal charges, integer single/double/triple bonds.
 * Do not silently infer topology or claim to preserve other CJSON scientific properties. */
export function exportCJSON(doc){
 const indices=new Map(doc.atoms.map((a,i)=>[a.id,i]));
 const payload={chemicalJson:1,name:String(doc.name||'MolecolaLab').slice(0,80),
  atoms:{elements:{number:doc.atoms.map(a=>element(a.symbol).Z)},
   coords:{'3d':doc.atoms.flatMap(a=>[a.x,a.y,a.z])},
   formalCharges:doc.atoms.map(a=>a.charge||0)},
  bonds:{connections:{index:doc.bonds.flatMap(b=>[indices.get(b.a),indices.get(b.b)])},
   order:doc.bonds.map(b=>b.order)}};
 return JSON.stringify(payload,null,2)+'\n';
}
export function importCJSON(text){
 let x;
 try{x=JSON.parse(text);}catch{throw Error('CJSON: JSON non valido');}
 if(!x||typeof x!=='object'||x.chemicalJson!==1||!x.atoms)
  throw Error('CJSON: richiesto chemicalJson: 1 e un oggetto atoms');
 if(x.unitCell||x.atoms.coords?.['3dFractional'])
  throw Error('CJSON: celle periodiche e coordinate frazionarie non supportate; il file originale resta invariato');
 const numbers=x.atoms.elements?.number,coords=x.atoms.coords?.['3d'];
 if(!Array.isArray(numbers)||!Array.isArray(coords))throw Error('CJSON: elementi o coordinate 3D mancanti');
 const n=numbers.length;
 if(n>6000||coords.length!==3*n)throw Error('CJSON: dimensione elenco coordinate incoerente o molecola troppo grande');
 if(!coords.every(v=>typeof v==='number'&&Number.isFinite(v)&&Math.abs(v)<=1000))throw Error('CJSON: coordinate cartesiane non valide (limite ±1000 Å)');
 const charges=x.atoms.formalCharges??Array(n).fill(0);
 if(!Array.isArray(charges)||charges.length!==n||charges.some(v=>!Number.isInteger(v)||Math.abs(v)>8))
  throw Error('CJSON: cariche formali non valide');
 const total=charges.reduce((sum,v)=>sum+v,0);
 if(x.properties?.totalCharge!==undefined && (!Number.isInteger(x.properties.totalCharge)||x.properties.totalCharge!==total))
  throw Error('CJSON: carica totale non coerente con le cariche formali; interrompo per evitare perdite');
 const d=blank();d.name=typeof x.name==='string'?(x.name.slice(0,80)||'Molecola CJSON'):'Molecola CJSON';
 for(let i=0;i<n;i++){
  const z=numbers[i];if(!Number.isInteger(z)||!atomicElement(z))throw Error('CJSON: numero atomico non valido');
  const id=addAtom(d,atomicElement(z).symbol,...coords.slice(3*i,3*i+3));
  findAtom(d,id).charge=charges[i];
 }
 if(x.bonds!==undefined){
  const links=x.bonds.connections?.index,orders=x.bonds.order;
  if(!Array.isArray(links)||!Array.isArray(orders)||orders.length>20000||links.length!==orders.length*2)
   throw Error('CJSON: indici dei legami non validi');
  const seen=new Set();
  for(let i=0;i<orders.length;i++){
   const a=links[2*i],b=links[2*i+1],order=orders[i];
   if(!Number.isInteger(a)||!Number.isInteger(b)||a<0||b<0||a>=n||b>=n||a===b||![1,2,3].includes(order))
    throw Error('CJSON: legame non supportato (aromatico/frazionario incluso) o indice non valido');
   const key=[Math.min(a,b),Math.max(a,b)].join('-');
   if(seen.has(key))throw Error('CJSON: legame duplicato');seen.add(key);
   setBond(d,a+1,b+1,order);
  }
 }
 return d;
}
/** CJSON can contain spectroscopy, orbitals, unit cells etc. Inform before discarding them. */
export function unhandledCJSONFields(text){
 const x=typeof text==='string'?JSON.parse(text):text;
 if(!x||typeof x!=='object'||x.chemicalJson!==1)return [];
 const extra=[];
 const topKnown=new Set(['chemicalJson','name','atoms','bonds']);
 for(const key of Object.keys(x))if(!topKnown.has(key)){
  if(key==='properties'&&x.properties&&Object.keys(x.properties).every(k=>k==='totalCharge'))continue;
  extra.push(key);
 }
 const atomKnown=new Set(['elements','coords','formalCharges']);
 if(x.atoms&&typeof x.atoms==='object')for(const key of Object.keys(x.atoms))if(!atomKnown.has(key))extra.push('atoms.'+key);
 const coordsKnown=new Set(['3d']);
 if(x.atoms?.coords)for(const key of Object.keys(x.atoms.coords))if(!coordsKnown.has(key))extra.push('atoms.coords.'+key);
 if(x.bonds&&typeof x.bonds==='object')for(const key of Object.keys(x.bonds))if(!['connections','order'].includes(key))extra.push('bonds.'+key);
 return extra;
}

/** Change a torsion A-B-C-D by rotating the entire D-side of the acyclic B-C single bond.
 * The B-C axis remains fixed; ring bonds and multiple bonds are intentionally rejected.
 * This is a coordinate transform, not a force-field optimization. */
export function setDihedral(doc,aId,bId,cId,dId,target){
 if(!Number.isFinite(target)||target< -180||target>180)throw Error('Diedro richiesto tra −180° e +180°');
 const ids=[aId,bId,cId,dId];
 if(new Set(ids).size!==4)throw Error('Seleziona quattro atomi distinti A–B–C–D');
 const atoms=ids.map(id=>findAtom(doc,id));
 if(atoms.some(a=>!a))throw Error('Uno o più atomi non esistono');
 if(!bondBetween(doc,aId,bId)||!bondBetween(doc,cId,dId))throw Error('Gli atomi devono formare una catena A–B–C–D');
 const pivot=bondBetween(doc,bId,cId);
 if(!pivot)throw Error('Il legame centrale B–C non esiste');
 if(pivot.order!==1)throw Error('Si possono ruotare solo legami centrali singoli');
 // Remove B–C from the graph. The C-side must be disconnected from B.
 const adjacency=new Map(doc.atoms.map(a=>[a.id,[]]));
 for(const bond of doc.bonds){
  if(bond.id===pivot.id)continue;
  adjacency.get(bond.a).push(bond.b);adjacency.get(bond.b).push(bond.a);
 }
 const moving=new Set([cId]),queue=[cId];
 for(let i=0;i<queue.length;i++)for(const other of adjacency.get(queue[i]))if(!moving.has(other)){moving.add(other);queue.push(other);}
 if(moving.has(bId))throw Error('Il legame centrale appartiene a un anello: rotazione bloccata');
 if(!moving.has(dId)||moving.has(aId))throw Error('La catena selezionata non definisce il lato rotante');
 const old=dihedral(...atoms);
 if(!Number.isFinite(old))throw Error('Diedro non definito: tre atomi collineari o atomi coincidenti');
 const [b,c]=[atoms[1],atoms[2]],axis=sub(c,b),norm=mag(axis);
 if(norm<1e-8)throw Error('Asse del legame centrale degenerato');
 const u=axis.map(v=>v/norm);
 const delta=((target-old+540)%360)-180;
 const t=-delta*Math.PI/180,co=Math.cos(t),si=Math.sin(t);
 const moved=[];
 for(const id of moving){if(id===cId)continue;
  const p=findAtom(doc,id),v=sub(p,b),pr=cross(u,v),project=dot(u,v);
  moved.push([p,{x:b.x+v[0]*co+pr[0]*si+u[0]*project*(1-co),
    y:b.y+v[1]*co+pr[1]*si+u[1]*project*(1-co),
    z:b.z+v[2]*co+pr[2]*si+u[2]*project*(1-co)}]);
 }
 // All positions are computed before mutating the document, so failures are atomic.
 for(const [p,q] of moved)setAtomPosition(doc,p.id,q.x,q.y,q.z);
 const achieved=dihedral(...ids.map(id=>findAtom(doc,id)));
 if(!Number.isFinite(achieved)||Math.abs((((achieved-target+540)%360)-180))>0.001)
  throw Error('Impossibile raggiungere il diedro richiesto');
 return {before:old,after:achieved,moved:moved.length};
}

export function exportProject(doc){return JSON.stringify({format:'molecolalab-v1',molecule:doc},null,2)+'\n';}
export function importProject(text){const x=JSON.parse(text);if(x.format!=='molecolalab-v1'||!x.molecule||!Array.isArray(x.molecule.atoms)||!Array.isArray(x.molecule.bonds))throw Error('File di progetto non valido');
 const d=x.molecule;if(d.atoms.length>6000||d.bonds.length>20000)throw Error('Progetto troppo grande');
 const ids=new Set();for(const a of d.atoms){if(!element(a.symbol)||!Number.isInteger(a.id)||a.id<1||ids.has(a.id)||![a.x,a.y,a.z].every(Number.isFinite)||!Number.isInteger(a.charge??0)||Math.abs(a.charge??0)>8)throw Error('Atomo nel progetto non valido');ids.add(a.id);}
 const bondIds=new Set(),pairs=new Set();
 for(const b of d.bonds){const pair=[b.a,b.b].sort((x,y)=>x-y).join(':');if(!ids.has(b.a)||!ids.has(b.b)||b.a===b.b||!Number.isInteger(b.id)||b.id<1||bondIds.has(b.id)||pairs.has(pair)||![1,2,3].includes(b.order))throw Error('Legame nel progetto non valido');bondIds.add(b.id);pairs.add(pair);}
 if(typeof d.name!=='string'||d.name.length>120)throw Error('Nome del progetto non valido');
 return clone(d);
}
function molecule(name,atomSpec,bondSpec){const d=blank();d.name=name;for(const [s,x,y,z] of atomSpec)addAtom(d,s,x,y,z);for(const [a,b,o=1] of bondSpec)setBond(d,a,b,o);return d;}
export const PRESETS={
 // E/Z examples: explicit hydrogens, generated offline from RDKit ETKDG + MMFF94.
 ezbutE:()=>molecule("(E)-but-2-ene",[["C",1.94965,-0.00545,-0.03846],["C",0.50963,-0.16136,-0.40525],["C",-0.50963,0.16136,0.40525],["C",-1.94965,0.00545,0.03846],["H",2.45178,-0.9768,-0.08005],["H",2.44326,0.66802,-0.74569],["H",2.07871,0.40287,0.96886],["H",0.30738,-0.56326,-1.39578],["H",-0.30737,0.56326,1.39578],["H",-2.45179,0.9768,0.08007],["H",-2.44325,-0.66803,0.74568],["H",-2.07871,-0.40285,-0.96887]],[[1,2,1],[2,3,2],[3,4,1],[1,5,1],[1,6,1],[1,7,1],[2,8,1],[3,9,1],[4,10,1],[4,11,1],[4,12,1]]),
 ezbutZ:()=>molecule("(Z)-but-2-ene",[["C",1.37166,-0.24791,-0.8639],["C",0.65607,-0.69438,0.3696],["C",-0.53046,-0.27125,0.8332],["C",-1.43747,0.75387,0.23368],["H",0.8328,0.51492,-1.42917],["H",2.34964,0.16263,-0.59369],["H",1.53358,-1.1045,-1.52575],["H",1.17899,-1.45612,0.94735],["H",-0.89646,-0.71598,1.75827],["H",-1.59398,1.56901,0.94715],["H",-1.05429,1.18785,-0.69186],["H",-2.41006,0.30187,0.01513]],[[1,2,1],[2,3,2],[3,4,1],[1,5,1],[1,6,1],[1,7,1],[2,8,1],[3,9,1],[4,10,1],[4,11,1],[4,12,1]]),
 ezdichloroE:()=>molecule("(E)-1,2-dicloroetene",[["Cl",2.10902,0.39232,0.0012],["C",0.56563,-0.34995,-0.00103],["C",-0.56563,0.34995,0.00103],["Cl",-2.10902,-0.39232,-0.00121],["H",0.62111,-1.43079,-0.00424],["H",-0.62111,1.43079,0.00424]],[[1,2,1],[2,3,2],[3,4,1],[2,5,1],[3,6,1]]),
 ezdichloroZ:()=>molecule("(Z)-1,2-dicloroetene",[["Cl",-1.62151,1.28497,0.00093],["C",-0.67065,-0.13927,0.00641],["C",0.66028,-0.18228,-0.00344],["Cl",1.70093,1.1776,-0.02366],["H",-1.28532,-1.03009,0.01914],["H",1.21628,-1.11094,0.00062]],[[1,2,1],[2,3,2],[3,4,1],[2,5,1],[3,6,1]]),
 ezacidE:()=>molecule("(E)-acido but-2-enoico",[["C",-2.05815,0.07675,-0.437],["C",-0.72314,-0.30301,0.11263],["C",0.43126,0.26085,-0.27366],["C",1.7059,-0.17027,0.31729],["O",1.8805,-1.02354,1.16514],["O",2.72016,0.52761,-0.21762],["H",-1.98557,0.8674,-1.19052],["H",-2.53191,-0.79401,-0.90036],["H",-2.70597,0.43257,0.3699],["H",-0.71673,-1.08616,0.86971],["H",0.4727,1.04227,-1.02252],["H",3.51095,0.16953,0.2365]],[[1,2,1],[2,3,2],[3,4,1],[4,5,2],[4,6,1],[1,7,1],[1,8,1],[1,9,1],[2,10,1],[3,11,1],[6,12,1]]),
 ezacidZ:()=>molecule("(Z)-acido but-2-enoico",[["C",-1.37214,-0.83226,-0.27197],["C",-0.93477,0.5964,-0.28295],["C",0.26656,1.10386,0.0382],["C",1.49606,0.39715,0.42802],["O",2.44345,0.93134,0.97764],["O",1.53026,-0.89356,0.07445],["H",-0.7899,-1.45485,0.4109],["H",-2.41255,-0.88874,0.06602],["H",-1.32193,-1.25235,-1.28059],["H",-1.71194,1.30468,-0.56705],["H",0.39617,2.18306,0.02918],["H",2.41073,-1.19474,0.37815]],[[1,2,1],[2,3,2],[3,4,1],[4,5,2],[4,6,1],[1,7,1],[1,8,1],[1,9,1],[2,10,1],[3,11,1],[6,12,1]]),
 brclfr:()=>molecule("Bromoclorofluorometano (R)", [["F",-0.0891,1.2879,0.5809],["C",-0.0591,0.1051,-0.0831],["Cl",-1.3038,-0.9705,0.5504],["Br",1.6876,-0.7115,0.0976],["H",-0.2356,0.2891,-1.1457]], [[1,2,1],[2,3,1],[2,4,1],[2,5,1]]),
 brclfs:()=>molecule("Bromoclorofluorometano (S)", [["F",-0.1601,1.1545,-0.8034],["C",-0.0556,0.1182,0.0661],["Cl",-1.3157,-1.0701,-0.2622],["Br",1.6931,-0.6998,-0.086],["H",-0.1618,0.4972,1.0855]], [[1,2,1],[2,3,1],[2,4,1],[2,5,1]]),
 butanolr:()=>molecule("(R)-butan-2-olo", [["C",-1.371,0.9148,0.0292],["C",-0.7139,-0.4581,0.0178],["O",-1.4344,-1.2706,-0.9088],["C",0.7579,-0.4225,-0.4028],["C",1.6485,0.3472,0.5586],["H",-2.4425,0.8222,0.2383],["H",-1.284,1.3989,-0.9498],["H",-0.9313,1.5704,0.7858],["H",-0.8153,-0.9169,1.0077],["H",-1.0443,-2.1609,-0.8834],["H",1.1299,-1.4517,-0.4845],["H",0.848,0.0018,-1.4106],["H",1.5583,-0.0447,1.5766],["H",2.696,0.258,0.2532],["H",1.398,1.4119,0.5728]], [[1,2,1],[2,3,1],[2,4,1],[4,5,1],[1,6,1],[1,7,1],[1,8,1],[2,9,1],[3,10,1],[4,11,1],[4,12,1],[5,13,1],[5,14,1],[5,15,1]]),
 butanols:()=>molecule("(S)-butan-2-olo", [["C",-1.366,0.982,-0.076],["C",-0.7474,-0.371,-0.3961],["O",-1.0809,-1.2736,0.6544],["C",0.7723,-0.3153,-0.5863],["C",1.544,0.1319,0.6486],["H",-2.4574,0.8951,-0.0315],["H",-1.1094,1.7276,-0.8344],["H",-1.0492,1.348,0.9061],["H",-1.2045,-0.7548,-1.3155],["H",-0.6994,-2.1392,0.4293],["H",1.1288,-1.3173,-0.8554],["H",1.0153,0.349,-1.4237],["H",1.3313,-0.5117,1.5077],["H",2.6204,0.0853,0.4533],["H",1.302,1.1638,0.9194]], [[1,2,1],[2,3,1],[2,4,1],[4,5,1],[1,6,1],[1,7,1],[1,8,1],[2,9,1],[3,10,1],[4,11,1],[4,12,1],[5,13,1],[5,14,1],[5,15,1]]),
 lactic_s:()=>molecule("(S)-acido lattico", [["C",-1.3873,-0.5477,0.4319],["C",-0.2772,0.4488,0.1501],["O",-0.5693,1.0994,-1.086],["C",1.0858,-0.245,0.1214],["O",1.4579,-1.1875,0.799],["O",1.9318,0.3471,-0.7563],["H",-1.2422,-1.0482,1.3943],["H",-1.4293,-1.3158,-0.3488],["H",-2.3617,-0.0478,0.4398],["H",-0.2439,1.2146,0.932],["H",0.2817,1.4663,-1.3974],["H",2.7536,-0.184,-0.6799]], [[1,2,1],[2,3,1],[2,4,1],[4,5,2],[4,6,1],[1,7,1],[1,8,1],[1,9,1],[2,10,1],[3,11,1],[6,12,1]]),
 alanine_s:()=>molecule("(S)-alanina", [["N",0.18,1.4649,-0.6188],["C",-0.1994,0.4399,0.3766],["C",-1.3228,-0.4412,-0.1528],["C",1.0,-0.4376,0.7426],["O",1.1749,-1.0126,1.8057],["O",1.875,-0.61,-0.2756],["H",-0.6538,1.8279,-1.0761],["H",0.7339,1.0032,-1.3451],["H",-0.5208,0.9603,1.285],["H",-1.0231,-0.969,-1.0656],["H",-1.6002,-1.2008,0.5865],["H",-2.2166,0.1504,-0.3781],["H",2.5729,-1.1753,0.1158]], [[1,2,1],[2,3,1],[2,4,1],[4,5,2],[4,6,1],[1,7,1],[1,8,1],[2,9,1],[3,10,1],[3,11,1],[3,12,1],[6,13,1]]),
 ditartaric:()=>molecule("(2R,3R)-acido tartarico", [["O",1.3258,1.5044,0.4309],["C",1.3125,0.6469,-0.4353],["O",1.9152,0.8692,-1.6298],["C",0.7183,-0.7615,-0.3588],["O",0.426,-1.274,-1.6678],["C",-0.5411,-0.8212,0.508],["O",-0.1698,-1.0054,1.8825],["C",-1.4319,0.4111,0.3352],["O",-1.6168,1.0681,-0.6748],["O",-2.091,0.7141,1.4813],["H",2.233,1.7959,-1.5606],["H",1.4993,-1.4053,0.0625],["H",1.1081,-0.8913,-2.2567],["H",-1.1521,-1.6879,0.2298],["H",-0.9321,-0.6819,2.4047],["H",-2.6036,1.5186,1.249]], [[1,2,2],[2,3,1],[2,4,1],[4,5,1],[4,6,1],[6,7,1],[6,8,1],[8,9,2],[8,10,1],[3,11,1],[4,12,1],[5,13,1],[6,14,1],[7,15,1],[10,16,1]]),
 water:()=>molecule('Acqua (H₂O)',[['O',0,0,0],['H',-.757,.586,0],['H',.757,.586,0]],[[1,2],[1,3]]),
 ethane:()=>molecule('Etano (C₂H₆) — geometria MMFF94 di riferimento',[['C',0.7445,-0.1272,0.0342],['C',-0.7445,0.1272,-0.0342],['H',1.2514,0.3243,-0.8238],['H',1.166,0.3017,0.9482],['H',0.9515,-1.2015,0.0305],['H',-1.166,-0.3017,-0.9482],['H',-1.2514,-0.3243,0.8238],['H',-0.9515,1.2015,-0.0305]],[[1,2],[1,3],[1,4],[1,5],[2,6],[2,7],[2,8]]),
 butane:()=>molecule('Butano (C₄H₁₀) — conformero MMFF94 di riferimento',[['C',1.3714,0.6264,-0.6246],['C',0.8012,-0.6208,0.033],['C',-0.4505,-0.3503,0.8682],['C',-1.6397,0.1166,0.043],['H',1.576,1.4038,0.1182],['H',0.6859,1.0341,-1.3733],['H',2.3116,0.3869,-1.1319],['H',0.5818,-1.3704,-0.7362],['H',1.5688,-1.0525,0.6864],['H',-0.7284,-1.2749,1.3882],['H',-0.2293,0.3941,1.6419],['H',-1.8678,-0.5959,-0.7559],['H',-2.526,0.207,0.6793],['H',-1.455,1.0959,-0.4082]],[[1,2],[2,3],[3,4],[1,5],[1,6],[1,7],[2,8],[2,9],[3,10],[3,11],[4,12],[4,13],[4,14]]),
  methane:()=>{const p=.629;return molecule('Metano (CH₄)',[['C',0,0,0],['H',p,p,p],['H',-p,-p,p],['H',-p,p,-p],['H',p,-p,-p]],[[1,2],[1,3],[1,4],[1,5]]);},
 ammonia:()=>molecule('Ammoniaca (NH₃)',[['N',0,0,.32],['H',.946,0,0],['H',-.473,.819,0],['H',-.473,-.819,0]],[[1,2],[1,3],[1,4]]),
 co2:()=>molecule('Diossido di carbonio (CO₂)',[['O',-1.16,0,0],['C',0,0,0],['O',1.16,0,0]],[[1,2,2],[2,3,2]]),
 ethene:()=>molecule('Etene (C₂H₄)',[['C',-.67,0,0],['C',.67,0,0],['H',-1.23,.93,0],['H',-1.23,-.93,0],['H',1.23,.93,0],['H',1.23,-.93,0]],[[1,2,2],[1,3],[1,4],[2,5],[2,6]]),
 ethyne:()=>molecule('Etino (C₂H₂)',[['C',-.60,0,0],['C',.60,0,0],['H',-1.66,0,0],['H',1.66,0,0]],[[1,2,3],[1,3],[2,4]]),
 formaldehyde:()=>molecule('Metanale (CH₂O)',[['C',0,0,0],['O',1.22,0,0],['H',-.57,.92,0],['H',-.57,-.92,0]],[[1,2,2],[1,3],[1,4]]),
 benzene:()=>{const a=[],b=[];for(let i=0;i<6;i++){const th=Math.PI*i/3;a.push(['C',1.39*Math.cos(th),1.39*Math.sin(th),0]);}for(let i=0;i<6;i++){const th=Math.PI*i/3;a.push(['H',2.48*Math.cos(th),2.48*Math.sin(th),0]);b.push([i+1,(i+1)%6+1,i%2===0?2:1],[i+1,i+7,1]);}return molecule('Benzene (C₆H₆) — una struttura di Kekulé',a,b);},
 bf3:()=>{const a=[['B',0,0,0]],b=[];for(let i=0;i<3;i++){const th=i*2*Math.PI/3;a.push(['F',1.3*Math.cos(th),1.3*Math.sin(th),0]);b.push([1,i+2]);}return molecule('Trifluoruro di boro (BF₃)',a,b);},
 sf6:()=>{const a=[['S',0,0,0],['F',1.56,0,0],['F',-1.56,0,0],['F',0,1.56,0],['F',0,-1.56,0],['F',0,0,1.56],['F',0,0,-1.56]];return molecule('Esafluoruro di zolfo (SF₆)',a,[[1,2],[1,3],[1,4],[1,5],[1,6],[1,7]]);},
 pcl5:()=>{const a=[['P',0,0,0]],b=[];for(let i=0;i<3;i++){const t=i*2*Math.PI/3;a.push(['Cl',2.02*Math.cos(t),2.02*Math.sin(t),0]);b.push([1,i+2]);}a.push(['Cl',0,0,2.14],['Cl',0,0,-2.14]);b.push([1,5],[1,6]);return molecule('Pentacloruro di fosforo (PCl₅) — geometria ideale',a,b);},
 hcn:()=>molecule('Acido cianidrico (HCN)',[['H',-2.22,0,0],['C',-1.16,0,0],['N',0,0,0]],[[1,2],[2,3,3]]),
 ammonium:()=>{const d=PRESETS.methane();d.name='Ione ammonio (NH₄⁺)';d.atoms[0].symbol='N';d.atoms[0].charge=1;return d;},
 carbonate:()=>{const a=[['C',0,0,0]];const b=[];for(let i=0;i<3;i++){const t=i*2*Math.PI/3;a.push(['O',1.29*Math.cos(t),1.29*Math.sin(t),0]);b.push([1,i+2,i===0?2:1]);}const d=molecule('Carbonato (CO₃²⁻) — una forma limite di risonanza',a,b);d.atoms[2].charge=-1;d.atoms[3].charge=-1;return d;},
};

import {
 ELEMENTS,element,blank,clone,findAtom,addAtom,setBond,eraseAtom,eraseBond,setAtomPosition,setAtomElement,setBondLength,estimateVSEPR,netCharge,setDihedral,
 dist,angle,dihedral,formula,centroid,center,suggestedPosition,hydrogenate,regularize,
 exportXYZ,importXYZ,exportMOL,importMOL,exportCJSON,importCJSON,unhandledCJSONFields,exportProject,importProject,PRESETS
} from './core.js';
import {analyzeCIP,summarizeCIP,mirrorMolecule,analyzeEZ,summarizeEZ,invertEZ} from './stereo.js';
const $=id=>document.getElementById(id);
const canvas=$('canvas'),ctx=canvas.getContext('2d',{alpha:true});
let doc=blank(), selected=null, selectedBond=null, pendingBond=null,measurement=[], tool='view',dirty=false;
let undo=[], redo=[], yaw=0.42, pitch=-0.17, zoom=1.0, mode='ballstick';
let pointerStart=null, dragging=false, dims={w:100,h:100},projected=[],stereoLabels=[],ezLabels=[];
const MAX_HISTORY=80;
const toolHelp={move:'Trascina un atomo per spostarlo nel piano dello schermo senza ricalcolare il resto della molecola. Per la profondità z usa le coordinate numeriche.',view:'Trascina la molecola per ruotarla, usa la rotellina per lo zoom. Clicca un atomo per selezionarlo.',atom:'Seleziona un elemento. Clicca sullo sfondo per aggiungerlo; prima clicca un atomo per creare un legame.',bond:'Clicca due atomi per inserire o modificare un legame dell’ordine scelto.',select:'Clicca un atomo per ispezionarlo; usa il pannello di destra per cambiare la carica formale.',measure:'Clicca 2, 3 o 4 atomi secondo la grandezza selezionata nel pannello destro.',delete:'Clicca un atomo o un legame per eliminarlo.'};
function toast(msg){const t=$('toast');t.textContent=msg;t.classList.add('show');clearTimeout(toast.timer);toast.timer=setTimeout(()=>t.classList.remove('show'),3700);}
function remember(before){if(JSON.stringify(before)===JSON.stringify(doc))return false;undo.push(before);if(undo.length>MAX_HISTORY)undo.shift();redo=[];dirty=true;return true;}
function changed(){render();refreshPanel();}
function applyMutator(fn){const before=clone(doc);try{fn();if(remember(before))changed();}catch(e){doc=before;changed();toast('Operazione non riuscita: '+e.message);}}
function clearSelection(){selected=null;selectedBond=null;pendingBond=null;measurement=[];changed();}
function replaceDoc(molecule){if(dirty&&doc.atoms.length&&!confirm('Sostituire la molecola? Le modifiche non esportate andranno perse.'))return false;doc=molecule;dirty=true;undo=[];redo=[];selected=null;selectedBond=null;pendingBond=null;measurement=[];resetCamera(false);changed();return true;}
function setTool(value){tool=value;pendingBond=null;measurement=[];
 document.querySelectorAll('[data-tool]').forEach(b=>{const active=b.dataset.tool===value;b.classList.toggle('active',active);b.setAttribute('aria-pressed',String(active));});
 canvas.classList.toggle('tool-mode',value!=='view');canvas.classList.toggle('move-mode',value==='move');$('toolHint').textContent=toolHelp[value];$('viewerMessage').textContent=toolHelp[value];changed();}
function initElements(){
 const s=$('elemSelect'),sel=$('selectedElement');for(const e of ELEMENTS){for(const target of [s,sel]){const opt=document.createElement('option');opt.value=e.symbol;opt.textContent=`${e.Z.toString().padStart(3,' ')} · ${e.symbol}`;target.appendChild(opt);}}s.value='C';
 for(const name of ['C','H','O','N','S','P','Cl','F','Br']){const b=document.createElement('button');b.textContent=name;b.dataset.symbol=name;b.addEventListener('click',()=>{s.value=name;refreshQuick();if(tool!=='atom')setTool('atom');});$('quickElements').appendChild(b);}
 s.addEventListener('change',refreshQuick);refreshQuick();
}
function refreshQuick(){document.querySelectorAll('[data-symbol]').forEach(b=>b.classList.toggle('active',b.dataset.symbol===$('elemSelect').value));}
function refreshPanel(){
 $('moleculeName').textContent=doc.name||'Nuova molecola';$('renameInput').value=doc.name||'';
 $('formula').textContent=formula(doc);
 const charge=netCharge(doc);$('chargeLabel').textContent=`Carica formale netta: ${charge>0?'+':''}${charge}`;
 $('counts').textContent=`${doc.atoms.length} atomi · ${doc.bonds.length} legami`;
 $('atomListCount').textContent=`(${doc.atoms.length})`;$('emptyState').hidden=doc.atoms.length>0;
 $('undoBtn').disabled=!undo.length;$('redoBtn').disabled=!redo.length;
 const a=selected===null?null:findAtom(doc,selected);
 const b=selectedBond===null?null:doc.bonds.find(t=>t.id===selectedBond);
 $('selectedTitle').textContent=a?`Atomo ${a.id} · ${a.symbol}`:b?`Legame ${b.id}`:'Nessuna selezione';
 $('atomControls').hidden=!a;$('bondControls').hidden=!b;
 const info=$('selectedInfo');info.replaceChildren();
 if(a){
  const n=doc.bonds.filter(t=>t.a===a.id||t.b===a.id).length;
  const entries=[['Numero atomico',element(a.symbol).Z],['Atomi direttamente legati',n],['Posizione x',a.x.toFixed(3)+' Å'],['Posizione y',a.y.toFixed(3)+' Å'],['Posizione z',a.z.toFixed(3)+' Å']];
  const dl=document.createElement('dl');for(const [k,v] of entries){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=k;dd.textContent=String(v);dl.append(dt,dd);}info.append(dl);
 }else if(b){
  const a1=findAtom(doc,b.a),a2=findAtom(doc,b.b);
  info.textContent=`${a1?.symbol||'?'}${b.a} — ${a2?.symbol||'?'}${b.b}: ${dist(a1,a2).toFixed(3)} Å`;
 }else info.textContent='Clicca un atomo o un legame nella vista 3D per modificarlo.';
 $('atomCharge').disabled=!a;$('atomCharge').value=String(a?.charge||0);
 $('selectedElement').disabled=!a;$('selectedElement').value=a?.symbol||'C';
 for(const axis of ['x','y','z']){const input=$('atom'+axis.toUpperCase());input.disabled=!a;input.value=a?Number(a[axis].toFixed(4)):'';}
 $('removeSelected').disabled=!a;
 $('selectedBondOrder').disabled=!b;$('selectedBondOrder').value=String(b?.order||1);
 $('selectedBondLength').disabled=!b;$('selectedBondLength').value=b?dist(findAtom(doc,b.a),findAtom(doc,b.b)).toFixed(3):'';
 $('applyBondLength').disabled=!b;$('removeSelectedBond').disabled=!b;
 const vsepr=$('vseprInfo');vsepr.replaceChildren();
 if(!a){vsepr.textContent='Seleziona un atomo centrale (almeno due vicini) per una stima Lewis/VSEPR.';}
 else {
  const prediction=estimateVSEPR(doc,a.id);
  if(!prediction?.supported){vsepr.textContent=prediction?.reason||'Previsione non disponibile.';}
  else {
   const value=document.createElement('strong');value.className='vsepr-value';value.textContent=`${prediction.formula} · ${prediction.molecularGeometry}`;
   const details=document.createElement('span');details.className='vsepr-detail';details.textContent=`Geometria elettronica: ${prediction.electronicGeometry}; ${prediction.electronDomains} domini, ${prediction.lonePairs} coppie solitarie stimate.`;
   const caveat=document.createElement('span');caveat.className='vsepr-detail';caveat.textContent=prediction.note;
   vsepr.append(value,details,caveat);
  }
 }
 refreshStereo(a);
 refreshEZ(b);
 const list=$('atomList');list.replaceChildren();if(!doc.atoms.length)list.textContent='Nessun atomo';
 for(const at of doc.atoms){
  const btn=document.createElement('button'),left=document.createElement('span'),right=document.createElement('span'),dot=document.createElement('span');
  dot.className='element-dot';dot.style.background=element(at.symbol).color;
  left.append(dot,document.createTextNode(`${at.id} · ${at.symbol}${at.charge?(at.charge>0?'+':'')+at.charge:''}`));
  right.textContent=`${at.x.toFixed(2)}, ${at.y.toFixed(2)}, ${at.z.toFixed(2)}`;
  btn.classList.toggle('active',at.id===selected);btn.append(left,right);btn.addEventListener('click',()=>{selected=at.id;selectedBond=null;changed();});list.append(btn);
 }
 $('measurementResult').replaceChildren();const n=Number($('measurementType').value);
 if(measurement.length>=n){const at=measurement.map(id=>findAtom(doc,id));if(at.every(Boolean)){
  if(n===4){['A','B','C','D'].forEach((letter,index)=>{$('torsion'+letter).value=String(measurement[index]);});
    const measured=dihedral(...at);if(Number.isFinite(measured))$('targetDihedral').value=measured.toFixed(2);
  }
  const v=n===2?dist(...at):n===3?angle(...at):dihedral(...at),unit=n===2?'Å':'°';
  const strong=document.createElement('strong');strong.textContent=Number.isFinite(v)?v.toFixed(n===2?3:2)+' '+unit:'Misura non definita';
  const desc=document.createElement('small');desc.textContent=`Atomi: ${measurement.join(' → ')}`;$('measurementResult').append(strong,desc);
 }}else $('measurementResult').textContent=`${measurement.length}/${n} atomi scelti. ${n===2?'Distanza':n===3?'Angolo di legame':'Angolo diedro'}.`;
}

/** R/S from topology and 3D coordinates, never from the view rotation. */
function refreshStereo(selectedAtom){
 const summary=$('stereoSummary'),inspector=$('stereoInfo');
 summary.replaceChildren();inspector.replaceChildren();
 if(!stereoLabels.length){summary.textContent='Nessun centro R/S assegnabile individuato.';}
 else{
  const assigned=stereoLabels.filter(x=>x.status==='assigned');
  const unresolved=stereoLabels.length-assigned.length;
  const header=document.createElement('small');header.className='muted';
  header.textContent=`${assigned.length} centri assegnati${unresolved?' · '+unresolved+' da verificare':''}`;summary.append(header);
  const box=document.createElement('div');box.className='stereo-badges';
  for(const x of stereoLabels){
   const b=document.createElement('button');b.type='button';b.className='stereo-chip '+(x.status==='assigned'?'assigned':'unresolved');
   b.textContent=x.status==='assigned'?`${x.symbol}${x.id}: ${x.descriptor}`:`${x.symbol}${x.id}: ?`;
   b.title=x.status==='assigned'?'Mostra le priorità CIP':x.reason;
   b.addEventListener('click',()=>{selected=x.id;selectedBond=null;changed();});box.append(b);
  }
  summary.append(box);
 }
 if(!selectedAtom){inspector.textContent='Seleziona un atomo C/Si/N⁺ tetraedrico. Clicca un centro R/S per leggere le priorità.';return;}
 const result=analyzeCIP(doc,selectedAtom.id);
 if(result.status==='assigned'){
  const heading=document.createElement('strong');heading.className='stereo-descriptor';heading.textContent=`${selectedAtom.symbol}${selectedAtom.id} · ${result.descriptor}`;
  inspector.append(heading);
  const title=document.createElement('div');title.className='stereo-priorities-title';title.textContent='Priorità CIP (1 = più alta)';inspector.append(title);
  const list=document.createElement('ol');list.className='stereo-priorities';
  for(const entry of result.priorities){const item=document.createElement('li');item.textContent=`${entry.symbol}${entry.atomId}`;list.append(item);}
  inspector.append(list);
  const note=document.createElement('small');note.className='stereo-note';note.textContent='La priorità dipende dalla connettività, R/S dalle coordinate tridimensionali. Ruotare la vista non cambia R/S.';inspector.append(note);
 }else{
  const first=document.createElement('strong');first.className='stereo-not-assigned';
  first.textContent=result.status==='not_stereogenic'?'Non stereogenico':'R/S non assegnato';inspector.append(first);
  const reason=document.createElement('span');reason.className='stereo-note';reason.textContent=result.reason;inspector.append(reason);
 }
}

/** E/Z inspector reflects actual 3D coordinates, not the camera orientation. */
function refreshEZ(selectedBondObject){
 const summary=$('ezSummary'),inspector=$('ezInfo');
 summary.replaceChildren();inspector.replaceChildren();
 const assigned=ezLabels.filter(x=>x.status==='assigned');
 const other=ezLabels.length-assigned.length;
 const heading=document.createElement('small');heading.className='muted';
 heading.textContent=`${assigned.length} doppi legami E/Z assegnati${other?' · '+other+' non assegnati':''}`;
 summary.append(heading);
 if(ezLabels.length){
  const chips=document.createElement('div');chips.className='stereo-badges';
  for(const entry of ezLabels){
   const button=document.createElement('button');button.type='button';
   button.className='stereo-chip '+(entry.status==='assigned'?'assigned':'unresolved');
   button.textContent=entry.status==='assigned'?`C${entry.firstId}=C${entry.secondId}: ${entry.descriptor}`:`C${entry.firstId}=C${entry.secondId}: ?`;
   button.title=entry.status==='assigned'?'Visualizza le priorità CIP':entry.reason;
   button.addEventListener('click',()=>{selectedBond=entry.bondId;selected=null;changed();});
   chips.append(button);
  }
  summary.append(chips);
 }else{
  const info=document.createElement('span');info.textContent='Nessun doppio legame C=C individuato.';summary.append(info);
 }
 let result=null;
 if(selectedBondObject)result=analyzeEZ(doc,selectedBondObject.id);
 const btn=$('invertEZBtn');
 btn.disabled=!selectedBondObject||result?.status!=='assigned';
 if(!selectedBondObject){inspector.textContent='Seleziona il centro di un doppio legame C=C o premi il badge della sua configurazione.';return;}
 if(result.status==='assigned'){
  const title=document.createElement('strong');title.className='stereo-descriptor';
  title.textContent=`C${result.firstId}=C${result.secondId} · ${result.descriptor}`;inspector.append(title);
  for(const end of result.ends){
   const sub=document.createElement('div');sub.className='stereo-priorities-title';
   sub.textContent=`Priorità CIP sul carbonio ${end.endId}`;inspector.append(sub);
   const list=document.createElement('ol');list.className='stereo-priorities ez-priorities';
   for(const lig of end.ranked){const item=document.createElement('li');item.textContent=`${lig.symbol}${lig.atomId}`;list.append(item);}
   inspector.append(list);
  }
  const note=document.createElement('small');note.className='stereo-note';
  note.textContent=result.descriptor==='E'?'E = i due sostituenti prioritari sono su lati opposti.':'Z = i due sostituenti prioritari sono dallo stesso lato.';
  inspector.append(note);
 }else{
  const title=document.createElement('strong');title.className='stereo-not-assigned';title.textContent='E/Z non assegnato';inspector.append(title);
  const note=document.createElement('small');note.className='stereo-note';note.textContent=result.reason;inspector.append(note);
 }
}

function resetCamera(refresh=true){
 yaw=.42;pitch=-.17;
 if(doc.atoms.length){const c=centroid(doc),radius=Math.max(...doc.atoms.map(a=>Math.hypot(a.x-c.x,a.y-c.y,a.z-c.z)));
  zoom=Math.min(3.1,Math.max(.55,3.4/Math.max(1,radius)));
 }else zoom=1;
 if(refresh)render();
}
function cameraPoint(atom){const c=centroid(doc),x=atom.x-c.x,y=atom.y-c.y,z=atom.z-c.z;
 const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);
 const xr=x*cy-z*sy,zr=x*sy+z*cy,yr=y*cp-zr*sp,zcam=y*sp+zr*cp;
 const scale=(Math.min(dims.w,dims.h)/9)*zoom*Math.max(.3,Math.min(3.5,10/(10-zcam)));
 return {id:atom.id,x:dims.w/2+xr*scale,y:dims.h/2-yr*scale,z:zcam,s:scale};
}
function unproject(sx,sy,zcam=0,referenceCenter=null){const c=referenceCenter||centroid(doc),cy=Math.cos(yaw),syw=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);
 const scale=(Math.min(dims.w,dims.h)/9)*zoom*Math.max(.3,Math.min(3.5,10/(10-zcam)));
 const xr=(sx-dims.w/2)/scale,yr=(dims.h/2-sy)/scale;
 const y=yr*cp+zcam*sp,zr=-yr*sp+zcam*cp;
 return {x:c.x+xr*cy+zr*syw,y:c.y+y,z:c.z-xr*syw+zr*cy};
}
function atomRadius(a,p){const r=element(a.symbol).radius;
 if(mode==='spacefill')return Math.max(8,Math.min(75,p.s*(.5+r*.34)));
 if(mode==='sticks')return Math.max(4,Math.min(14,p.s*.12));
 if(mode==='wire')return Math.max(3,Math.min(9,p.s*.085));
 return Math.max(8,Math.min(32,p.s*(.12+r*.12)));
}
const rgb=(hex)=>hex.match(/[a-f\d]{2}/gi).map(s=>parseInt(s,16));
function shade(hex,f){const [r,g,b]=rgb(hex);return `rgb(${Math.min(255,Math.max(0,r*f))|0},${Math.min(255,Math.max(0,g*f))|0},${Math.min(255,Math.max(0,b*f))|0})`;}
function drawBond(b){const a=projected.find(p=>p.id===b.a),c=projected.find(p=>p.id===b.b);if(!a||!c)return;
 const dx=c.x-a.x,dy=c.y-a.y,l=Math.hypot(dx,dy)||1,nx=-dy/l,ny=dx/l;
 const width=mode==='wire'?1.4:mode==='sticks'?8:mode==='spacefill'?3.5:5.5;
 const shift=b.order===1?[0]:b.order===2?[-.85,.85]:[-1.5,0,1.5];
 const aa=findAtom(doc,b.a),cc=findAtom(doc,b.b),ca=element(aa.symbol).color,cb=element(cc.symbol).color;
 if(selectedBond===b.id){ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(c.x,c.y);ctx.strokeStyle='#78ebdb';ctx.lineWidth=width+9;ctx.lineCap='round';ctx.stroke();}
 ctx.lineCap='round';for(const k of shift){const d=k*(mode==='wire'?3:width*.85);
  ctx.beginPath();ctx.moveTo(a.x+nx*d,a.y+ny*d);ctx.lineTo((a.x+c.x)/2+nx*d,(a.y+c.y)/2+ny*d);ctx.strokeStyle=mode==='wire'?shade(ca,.86):shade(ca,.9);ctx.lineWidth=width;ctx.stroke();
  ctx.beginPath();ctx.moveTo((a.x+c.x)/2+nx*d,(a.y+c.y)/2+ny*d);ctx.lineTo(c.x+nx*d,c.y+ny*d);ctx.strokeStyle=mode==='wire'?shade(cb,.86):shade(cb,.9);ctx.lineWidth=width;ctx.stroke();
 }
}
function drawAtom(a,p){const r=atomRadius(a,p),color=element(a.symbol).color;
 if(selected===a.id||measurement.includes(a.id)||pendingBond===a.id){ctx.beginPath();ctx.arc(p.x,p.y,r+6,0,Math.PI*2);ctx.strokeStyle='#78ebdb';ctx.lineWidth=2.4;ctx.stroke();}
 const gradient=ctx.createRadialGradient(p.x-r*.33,p.y-r*.42,Math.max(1,r*.02),p.x+r*.14,p.y+r*.14,r*1.26);
 gradient.addColorStop(0,shade(color,1.55));gradient.addColorStop(.36,shade(color,1.13));gradient.addColorStop(.75,shade(color,.8));gradient.addColorStop(1,shade(color,.5));
 ctx.fillStyle=gradient;ctx.beginPath();ctx.arc(p.x,p.y,r,0,Math.PI*2);ctx.fill();
 if(mode!=='wire'&&r>=9){ctx.fillStyle=a.symbol==='H'||a.symbol==='He'?'#183047':'#f7fbff';ctx.font=`${Math.max(9,Math.min(16,r*.72))|0}px system-ui`;ctx.textAlign='center';ctx.textBaseline='middle';ctx.fillText(a.symbol,p.x,p.y+.5);}
}
function render(){if(!ctx)return;const rect=canvas.getBoundingClientRect();const w=Math.max(1,rect.width),h=Math.max(1,rect.height);
 const pixelRatio=Math.min(window.devicePixelRatio||1,2);if(canvas.width!==Math.round(w*pixelRatio)||canvas.height!==Math.round(h*pixelRatio)){canvas.width=Math.round(w*pixelRatio);canvas.height=Math.round(h*pixelRatio);}dims={w,h};ctx.setTransform(pixelRatio,0,0,pixelRatio,0,0);ctx.clearRect(0,0,w,h);
 projected=doc.atoms.map(cameraPoint);
 const items=[...doc.bonds.map(b=>({type:'bond',depth:((projected.find(p=>p.id===b.a)?.z||0)+(projected.find(p=>p.id===b.b)?.z||0))/2,obj:b})),...doc.atoms.map(a=>({type:'atom',depth:projected.find(p=>p.id===a.id)?.z||0,obj:a}))];
 items.sort((a,b)=>a.depth-b.depth|| (a.type==='bond'?-1:1));for(const entry of items)if(entry.type==='bond')drawBond(entry.obj);else drawAtom(entry.obj,projected.find(p=>p.id===entry.obj.id));
 stereoLabels=summarizeCIP(doc);
 ezLabels=summarizeEZ(doc);
 if($('showRSLabels').checked){
  ctx.textAlign='center';ctx.textBaseline='middle';
  for(const center of stereoLabels.filter(c=>c.status==='assigned')){
   const p=projected.find(pt=>pt.id===center.id);if(!p)continue;
   const at=findAtom(doc,center.id);const offset=atomRadius(at,p)+12;
   const x=p.x+offset,y=p.y-offset;
   ctx.fillStyle=center.descriptor==='R'?'#e6b65b':'#73dcdb';ctx.beginPath();ctx.arc(x,y,13,0,2*Math.PI);ctx.fill();
   ctx.fillStyle='#101d30';ctx.font='800 13px system-ui';ctx.fillText(center.descriptor,x,y+.4);
  }
  if(selected!=null){
   const active=stereoLabels.find(item=>item.id===selected&&item.status==='assigned');
   if(active)for(const priority of active.priorities){
    const p=projected.find(v=>v.id===priority.atomId);if(!p)continue;
    const a=findAtom(doc,priority.atomId);const off=atomRadius(a,p)+10;
    const x=p.x+off*.66,y=p.y+off*.65;
    ctx.beginPath();ctx.arc(x,y,10,0,2*Math.PI);ctx.fillStyle='#fff1cf';ctx.fill();
    ctx.fillStyle='#0f1c2b';ctx.font='800 11px system-ui';ctx.fillText(priority.rank,x,y);
   }
  }
 }
 if($('showEZLabels').checked){
  ctx.textAlign='center';ctx.textBaseline='middle';
  for(const edge of ezLabels.filter(item=>item.status==='assigned')){
   const p=projected.find(v=>v.id===edge.firstId),q=projected.find(v=>v.id===edge.secondId);
   if(!p||!q)continue;
   const x=(p.x+q.x)/2,y=(p.y+q.y)/2;
   const vx=q.x-p.x,vy=q.y-p.y,len=Math.hypot(vx,vy)||1;
   const cx=x-vy/len*23,cy=y+vx/len*23;
   ctx.fillStyle=edge.descriptor==='E'?'#f4c687':'#b3a3ff';ctx.beginPath();ctx.arc(cx,cy,12,0,2*Math.PI);ctx.fill();
   ctx.fillStyle='#152132';ctx.font='800 13px system-ui';ctx.fillText(edge.descriptor,cx,cy);
   if(edge.bondId===selectedBond){
    for(const end of edge.ends){
     for(const rank of end.ranked){
      const pos=projected.find(v=>v.id===rank.atomId);if(!pos)continue;
      const atom=findAtom(doc,rank.atomId),offset=atomRadius(atom,pos)+11;
      const xx=pos.x+offset*.7,yy=pos.y+offset*.7;
      ctx.beginPath();ctx.arc(xx,yy,10,0,2*Math.PI);ctx.fillStyle='#fdf3ca';ctx.fill();
      ctx.fillStyle='#101d30';ctx.font='800 11px system-ui';ctx.fillText(rank.rank,xx,yy);
     }
    }
   }
  }
 }

}
function localPos(event){const r=canvas.getBoundingClientRect();return {x:event.clientX-r.left,y:event.clientY-r.top};}
function pick(x,y){let id=null,best=-Infinity;for(const p of projected){const a=findAtom(doc,p.id),radius=atomRadius(a,p)+7,d=Math.hypot(x-p.x,y-p.y);if(d<=Math.max(12,radius)){const score=-d/Math.max(10,radius)+p.z*.04;if(score>best){best=score;id=p.id;}}}
 if(id!==null)return {type:'atom',id};
 for(const b of doc.bonds){const a=projected.find(p=>p.id===b.a),c=projected.find(p=>p.id===b.b);if(!a||!c)continue;const vx=c.x-a.x,vy=c.y-a.y,lengthSq=vx*vx+vy*vy;if(!lengthSq)continue;
 const t=Math.max(0,Math.min(1,((x-a.x)*vx+(y-a.y)*vy)/lengthSq)),d=Math.hypot(x-(a.x+t*vx),y-(a.y+t*vy));if(d<9)return {type:'bond',id:b.id};
 }
 return null;
}
function clickScene(x,y){const p=pick(x,y);if(tool==='view'||tool==='select'||tool==='move'){selected=p?.type==='atom'?p.id:null;selectedBond=p?.type==='bond'?p.id:null;changed();return;}
 if(tool==='atom'){if(p?.type==='atom'){selected=p.id;selectedBond=null;toast(`Atomo ${p.id} scelto. Clicca uno spazio libero per collegarne un altro.`);changed();return;}
  let pos;if(selected&&findAtom(doc,selected)){
   const anchor=findAtom(doc,selected),cam=projected.find(p=>p.id===anchor.id),approx=unproject(x,y,cam?.z||0);
   let vec=[approx.x-anchor.x,approx.y-anchor.y,approx.z-anchor.z],m=Math.hypot(...vec);
   if(m<.3){pos=suggestedPosition(doc,anchor.id);}else{const len=element(anchor.symbol).radius+element($('elemSelect').value).radius;pos={x:anchor.x+vec[0]*len/m,y:anchor.y+vec[1]*len/m,z:anchor.z+vec[2]*len/m};}
  }else pos=unproject(x,y,0);
  let created=null;applyMutator(()=>{created=addAtom(doc,$('elemSelect').value,pos.x,pos.y,pos.z,selected&&findAtom(doc,selected)?selected:null,Number($('bondOrder').value));});if(created){selected=created;changed();}return;
 }
 if(tool==='bond'){if(p?.type!=='atom'){toast('Seleziona un atomo.');return;}
  if(pendingBond===null){pendingBond=p.id;selected=p.id;selectedBond=null;changed();toast('Seleziona il secondo atomo del legame.');return;}
  if(pendingBond!==p.id){const a=pendingBond;applyMutator(()=>setBond(doc,a,p.id,Number($('bondOrder').value)));selected=p.id;}
  pendingBond=null;changed();return;
 }
 if(tool==='measure'){if(p?.type!=='atom')return;
  const n=Number($('measurementType').value);if(measurement.length>=n)measurement=[];measurement.push(p.id);selected=p.id;selectedBond=null;changed();return;
 }
 if(tool==='delete'){if(!p)return;
  if(p.type==='atom')applyMutator(()=>eraseAtom(doc,p.id));else applyMutator(()=>eraseBond(doc,p.id));selected=null;selectedBond=null;changed();}
}
function setupPointer(){
 const pointers=new Map();let pinch=null,multiTouch=false;
 const pointerDistance=()=>{const p=[...pointers.values()];return Math.hypot(p[0].x-p[1].x,p[0].y-p[1].y);};
 canvas.addEventListener('pointerdown',e=>{
  if(e.pointerType==='mouse'&&e.button!==0&&e.button!==1)return;
  canvas.setPointerCapture(e.pointerId);const pos=localPos(e);pointers.set(e.pointerId,pos);
  if(pointers.size===1){
   pointerStart={...pos,yaw,pitch,button:e.button};dragging=false;
   if(tool==='move'&&e.button===0){const hit=pick(pos.x,pos.y);
    if(hit?.type==='atom'){
     const projectedAtom=projected.find(p=>p.id===hit.id);
     pointerStart.moving={id:hit.id,depth:projectedAtom.z,refCenter:centroid(doc),offsetX:pos.x-projectedAtom.x,offsetY:pos.y-projectedAtom.y,original:clone(doc)};
     selected=hit.id;selectedBond=null;changed();
    }
   }
  }else if(pointers.size>=2){
   if(pointerStart?.moving){doc=pointerStart.moving.original;changed();}
   multiTouch=true;pointerStart=null;pinch={distance:Math.max(5,pointerDistance()),zoom};
  }
 });
 canvas.addEventListener('pointermove',e=>{
  if(!pointers.has(e.pointerId))return;const pos=localPos(e);pointers.set(e.pointerId,pos);
  if(pointers.size>=2&&pinch){zoom=Math.min(4.5,Math.max(.25,pinch.zoom*pointerDistance()/pinch.distance));render();return;}
  if(multiTouch||!pointerStart)return;
  const dx=pos.x-pointerStart.x,dy=pos.y-pointerStart.y;
  if(Math.hypot(dx,dy)>5)dragging=true;
  if(dragging&&pointerStart.moving){
   const m=pointerStart.moving,point=unproject(pos.x-m.offsetX,pos.y-m.offsetY,m.depth,m.refCenter);
   try{setAtomPosition(doc,m.id,point.x,point.y,point.z);render();refreshPanel();}catch(err){toast(err.message);}
  }else if(dragging&&(tool==='view'||pointerStart.button===1)){
   yaw=pointerStart.yaw+dx*.008;pitch=Math.max(-Math.PI/2+.01,Math.min(Math.PI/2-.01,pointerStart.pitch+dy*.008));render();
  }
 });
 function endPointer(e,canceled=false){
  if(!pointers.has(e.pointerId))return;
  const wasMulti=multiTouch,wasDrag=dragging,pos=localPos(e),start=pointerStart;
  const didRotate=(tool==='view'||start?.button===1);
  pointers.delete(e.pointerId);pointerStart=null;dragging=false;
  if(start?.moving){
   if(canceled){doc=start.moving.original;changed();}
   else if(wasDrag){remember(start.moving.original);changed();}
   else if(!wasMulti){selected=start.moving.id;selectedBond=null;changed();}
  }else if(!wasMulti&&!canceled&&!wasDrag&&e.button===0)clickScene(pos.x,pos.y);
  else if(!wasMulti&&!canceled&&wasDrag&&!didRotate)toast('Trascina con «Sposta» oppure scegli «Ruota» per orientare la scena.');
  if(pointers.size===0){pinch=null;multiTouch=false;}
 }
 canvas.addEventListener('pointerup',e=>endPointer(e));
 canvas.addEventListener('pointercancel',e=>endPointer(e,true));
 canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.min(4.5,Math.max(.25,zoom*Math.exp(-e.deltaY*.001)));render();},{passive:false});
 canvas.addEventListener('dblclick',()=>resetCamera());
 new ResizeObserver(render).observe(canvas);
}

function exportFile(){let content,ext,mime='text/plain;charset=utf-8';const format=$('exportFormat').value;
 try{if(format==='mol'){content=exportMOL(doc);ext='mol';}else if(format==='xyz'){content=exportXYZ(doc);ext='xyz';}else if(format==='cjson'){content=exportCJSON(doc);ext='cjson';mime='application/json';}else{content=exportProject(doc);ext='json';mime='application/json';}}catch(e){toast(e.message);return;}
 const blob=new Blob([content],{type:mime}),url=URL.createObjectURL(blob),link=document.createElement('a');
 link.href=url;link.download=`molecolalab-${(doc.name||'molecola').normalize('NFD').replace(/[\u0300-\u036f]/g,'').toLowerCase().replace(/[^a-z0-9-]/g,'-').slice(0,36)}.${ext}`;
 document.body.appendChild(link);link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),1500);
 if(format==='mlab')dirty=false;toast(`Esportato ${ext.toUpperCase()}. I dati sono rimasti nel tuo browser.`);
}
function openFile(){const input=$('fileInput');input.value='';input.click();}
function guessLoader(text,name){const low=name.toLowerCase();if(low.endsWith('.cjson'))return importCJSON(text);if(low.endsWith('.json')||low.endsWith('.mlab')){const parsed=JSON.parse(text);return parsed?.chemicalJson===1?importCJSON(text):importProject(text);}
 if(low.endsWith('.xyz'))return importXYZ(text);
 if(low.endsWith('.mol')||low.endsWith('.sdf'))return importMOL(text.split(/\$\$\$\$/)[0]);
 throw Error('Formato file non riconosciuto: usa .mol, .sdf, .xyz, .json o .cjson.');}
function setup(){initElements();setupPointer();setTool('view');
 document.querySelectorAll('[data-tool]').forEach(b=>b.addEventListener('click',()=>setTool(b.dataset.tool)));
 $('presets').addEventListener('change',e=>{if(e.target.value&&PRESETS[e.target.value]){if(replaceDoc(PRESETS[e.target.value]()))toast('Modello di esempio caricato.');}e.target.value='';});
 $('startExample').addEventListener('click',()=>replaceDoc(PRESETS.water()));
 $('displayMode').addEventListener('change',e=>{mode=e.target.value;render();});
 $('newBtn').addEventListener('click',()=>replaceDoc(blank()));$('openBtn').addEventListener('click',openFile);$('saveBtn').addEventListener('click',exportFile);
 $('fileInput').addEventListener('change',async e=>{const file=e.target.files?.[0];if(!file)return;if(file.size>2e6){toast('File troppo grande (massimo 2 MB).');return;}
  try{const raw=await file.text();const model=guessLoader(raw,file.name);
   if(file.name.toLowerCase().endsWith('.cjson')||(file.name.toLowerCase().endsWith('.json')&&JSON.parse(raw)?.chemicalJson===1)){
    const extra=unhandledCJSONFields(raw);
    if(extra.length&&!confirm('Il file CJSON contiene dati non supportati ('+extra.slice(0,8).join(', ')+'). Verranno importati SOLO atomi, coordinate, cariche e legami; gli altri dati non saranno conservati. Conserva il file originale. Continuare?'))return;
   }
   if(replaceDoc(model))toast('File importato. Controlla cariche, connettività e proprietà scientifiche.');
  }catch(err){toast('Importazione fallita: '+err.message);}
 });
 $('undoBtn').addEventListener('click',undoStep);$('redoBtn').addEventListener('click',redoStep);
 $('addHsBtn').addEventListener('click',()=>{let n=0;applyMutator(()=>{n=hydrogenate(doc);});toast(n?`Aggiunti ${n} H secondo le valenze neutre semplificate.`:'Nessun H da aggiungere secondo la regola semplificata.');});
 $('regularizeBtn').addEventListener('click',()=>{if(doc.atoms.length>220){toast('Riordino grafico disponibile fino a 220 atomi in questa versione.');return;}applyMutator(()=>regularize(doc));toast('Riordinata la spaziatura (non è un’ottimizzazione energetica).');});
 $('centerBtn').addEventListener('click',()=>{applyMutator(()=>center(doc));toast('Coordinate centrate.');});
 $('resetCameraBtn').addEventListener('click',()=>resetCamera());
 $('showRSLabels').addEventListener('change',()=>render());
 $('showEZLabels').addEventListener('change',()=>render());
 $('invertEZBtn').addEventListener('click',()=>{if(selectedBond==null)return;let result=null;
  applyMutator(()=>{result=invertEZ(doc,selectedBond);});
  if(result)toast(`Creato isomero ${result.after} dal precedente ${result.before}. La modifica non rappresenta una rotazione fisica del doppio legame.`);
 });
 $('mirrorBtn').addEventListener('click',()=>{applyMutator(()=>mirrorMolecule(doc));toast('Molecola speculare creata: le configurazioni dei centri definiti vengono invertite.');});
 $('renameBtn').addEventListener('click',()=>{const name=$('renameInput').value.trim();if(!name){toast('Inserisci un nome non vuoto.');return;}applyMutator(()=>{doc.name=name.slice(0,80);});});
 $('clearSelection').addEventListener('click',clearSelection);
 $('selectedElement').addEventListener('change',e=>{if(selected!=null){applyMutator(()=>setAtomElement(doc,selected,e.target.value));toast('Elemento modificato: verifica legami, cariche e valenze.');}});
 $('applyCoords').addEventListener('click',()=>{if(selected==null)return;const raw=['X','Y','Z'].map(axis=>$('atom'+axis).value.trim());if(raw.some(v=>v==='')){toast('Inserisci tutte le coordinate.');return;}applyMutator(()=>setAtomPosition(doc,selected,...raw.map(Number)));});
 $('selectedBondOrder').addEventListener('change',e=>{const b=doc.bonds.find(b=>b.id===selectedBond);if(b)applyMutator(()=>setBond(doc,b.a,b.b,Number(e.target.value)));});
 $('applyBondLength').addEventListener('click',()=>{if(selectedBond==null)return;const raw=$('selectedBondLength').value.trim();if(!raw){toast('Inserisci una lunghezza in Å.');return;}applyMutator(()=>setBondLength(doc,selectedBond,Number(raw)));});
 $('removeSelectedBond').addEventListener('click',()=>{if(selectedBond==null)return;applyMutator(()=>eraseBond(doc,selectedBond));selectedBond=null;changed();});
 $('removeSelected').addEventListener('click',()=>{if(selected==null)return;applyMutator(()=>eraseAtom(doc,selected));selected=null;changed();});
 $('atomCharge').addEventListener('change',e=>{if(selected==null)return;applyMutator(()=>{const a=findAtom(doc,selected);if(a)a.charge=Number(e.target.value);});});
 $('readDihedralBtn').addEventListener('click',()=>{
  try{const ids=['A','B','C','D'].map(letter=>Number($('torsion'+letter).value.trim()));
    if(ids.some(id=>!Number.isSafeInteger(id)||id<=0))throw Error('Inserisci quattro numeri di atomo validi.');
    const atoms=ids.map(id=>findAtom(doc,id));if(atoms.some(a=>!a))throw Error('Uno o più atomi non esistono.');
    const value=dihedral(...atoms);
    if(!Number.isFinite(value))throw Error('Angolo diedro non definito per questa geometria.');
    setDihedral(clone(doc),...ids,value); // validate actual bond sequence, cycles and bond order, without mutation
    $('targetDihedral').value=value.toFixed(2);
    $('dihedralFeedback').textContent=`Diedro ${ids.join('–')}: ${value.toFixed(2)}°. Puoi ora modificarlo.`;
  }catch(e){$('dihedralFeedback').textContent=e.message;toast(e.message);}
 });
 $('applyDihedralBtn').addEventListener('click',()=>{
  const raw=$('targetDihedral').value.trim();if(!raw){toast('Indica l’angolo diedro richiesto.');return;}
  const ids=['A','B','C','D'].map(letter=>Number($('torsion'+letter).value.trim()));
  let result;
  applyMutator(()=>{result=setDihedral(doc,...ids,Number(raw));});
  if(result){$('dihedralFeedback').textContent=`Diedro ${ids.join('–')}: ${result.before.toFixed(2)}° → ${result.after.toFixed(2)}°. ${result.moved} atomi del frammento spostati.`;}
 });
 $('clearMeasurements').addEventListener('click',()=>{measurement=[];changed();});
 $('measurementType').addEventListener('change',()=>{measurement=[];changed();});
 $('helpBtn').addEventListener('click',()=>{
  const content=$('helpDialogContent');content.replaceChildren();
  const copy=$('instructionsBody').cloneNode(true);copy.removeAttribute('id');
  // Expand the first section for immediate orientation; the user can open the rest.
  content.append(copy);$('helpDialog').showModal();
 });
 $('closeHelp').addEventListener('click',()=>$('helpDialog').close());
 document.addEventListener('keydown',e=>{if(['INPUT','SELECT','TEXTAREA'].includes(document.activeElement?.tagName)||$('helpDialog').open)return;
  if((e.ctrlKey||e.metaKey)&&e.key.toLowerCase()==='z'){e.preventDefault();e.shiftKey?redoStep():undoStep();}
  else if((e.ctrlKey||e.metaKey)&&e.key.toLowerCase()==='y'){e.preventDefault();redoStep();}
  else if(e.key==='Escape'){e.preventDefault();clearSelection();}
  else if((e.key==='Delete'||e.key==='Backspace')&&(selected!=null||selectedBond!=null)){
   e.preventDefault();if(selected!=null)applyMutator(()=>eraseAtom(doc,selected));else applyMutator(()=>eraseBond(doc,selectedBond));clearSelection();
  }
 });
 window.addEventListener('beforeunload',e=>{if(dirty&&doc.atoms.length){e.preventDefault();e.returnValue='';}});
 changed();
}
function undoStep(){if(!undo.length)return;redo.push(clone(doc));doc=undo.pop();selected=null;selectedBond=null;pendingBond=null;measurement=[];dirty=true;changed();}
function redoStep(){if(!redo.length)return;undo.push(clone(doc));doc=redo.pop();selected=null;selectedBond=null;pendingBond=null;measurement=[];dirty=true;changed();}
setup();

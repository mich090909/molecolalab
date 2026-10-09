import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {PRESETS,dist,clone,bondBetween,exportProject,importProject,exportCJSON,importCJSON} from '../src/core.js';
import {analyzeEZ,summarizeEZ,invertEZ,analyzeCIP,mirrorMolecule} from '../src/stereo.js';
const models=JSON.parse(readFileSync(new URL('./fixtures/ez_rdkit_cases.json',import.meta.url),'utf8'));
const allDistances=doc=>doc.bonds.map(b=>dist(doc.atoms.find(a=>a.id===b.a),doc.atoms.find(a=>a.id===b.b)));
const distStable=(before,after)=>{const x=allDistances(before),y=allDistances(after);assert.equal(x.length,y.length);x.forEach((n,i)=>assert.ok(Math.abs(n-y[i])<1e-7,`bond ${i} distorted: ${n} -> ${y[i]}`));};
for(const [i,model] of models.entries()){
 const index=Object.keys(model.expected).map(Number)[0],expected=model.expected[index];
 test(`RDKit E/Z reference ${i+1} (${model.inputSmiles}) = ${expected}`,()=>{
  const result=analyzeEZ(model,index);
  assert.equal(result.status,'assigned',result.reason);
  assert.equal(result.descriptor,expected);
  assert.equal(result.ends.length,2);
  for(const side of result.ends){assert.deepEqual(side.ranked.map(k=>k.rank),[1,2]);assert.equal(new Set(side.ranked.map(k=>k.atomId)).size,2);}
 });
 test(`Inversion of E/Z reference ${i+1}: ${expected} -> ${expected==='E'?'Z':'E'} and back; lengths preserved`,()=>{
  const working=clone(model),first=clone(working),origBondCount=working.bonds.length;
  const result=invertEZ(working,index);
  assert.equal(result.before,expected);assert.equal(result.after,expected==='E'?'Z':'E');
  assert.equal(analyzeEZ(working,index).descriptor,result.after);
  assert.equal(working.bonds.length,origBondCount);
  distStable(first,working);
  const second=invertEZ(working,index);
  assert.equal(second.after,expected);
  assert.equal(analyzeEZ(working,index).descriptor,expected);
  distStable(first,working);
  for(let n=0;n<first.atoms.length;n++)for(const key of ['x','y','z'])
   assert.ok(Math.abs(first.atoms[n][key]-working.atoms[n][key])<1e-6,'double inversion changed atomic coordinates');
 });
 test(`Export/import preserves reference ${i+1} E/Z descriptor`,()=>{
  for(const saved of [exportProject(model),exportCJSON(model)]){
   const restored=saved.includes('molecolalab-v1')?importProject(saved):importCJSON(saved);
   assert.equal(analyzeEZ(restored,index).descriptor,expected);
  }
 });
}

test('Presets include three E/Z pairs with chemically distinct ligands',()=>{
 for(const pair of [['ezbutE','ezbutZ'],['ezdichloroE','ezdichloroZ'],['ezacidE','ezacidZ']]){
  const left=summarizeEZ(PRESETS[pair[0]]()),right=summarizeEZ(PRESETS[pair[1]]());
  assert.equal(left.length,1);assert.equal(right.length,1);
  assert.equal(left[0].descriptor,'E');assert.equal(right[0].descriptor,'Z');
 }
});
test('Ethene cannot have an E/Z descriptor because ligands coincide',()=>{
 const m=PRESETS.ethene(),bond=m.bonds.find(b=>b.order===2),r=analyzeEZ(m,bond.id);
 assert.equal(r.status,'not_stereogenic');
 assert.throws(()=>invertEZ(m,bond.id),/non assegnabile/);
});
test('Non-C=C double bond is refused, with no unsupported assignment',()=>{
 const m=PRESETS.formaldehyde();
 assert.equal(analyzeEZ(m,m.bonds.find(b=>b.order===2).id).status,'unsupported');
});
test('Triple and single bonds cannot be tagged E/Z',()=>{
 const m=PRESETS.ethyne();assert.equal(analyzeEZ(m,m.bonds[0].id).status,'not_double');
 assert.equal(analyzeEZ(m,m.bonds[1].id).status,'not_double');
});
test('Missing explicit H does not silently produce an E/Z descriptor',()=>{
 const m=PRESETS.ezbutE();m.atoms=m.atoms.filter(a=>a.symbol!=='H');m.bonds=m.bonds.filter(b=>m.atoms.some(a=>a.id===b.a)&&m.atoms.some(a=>a.id===b.b));
 assert.equal(analyzeEZ(m,m.bonds.find(b=>b.order===2).id).status,'unsupported');
});
test('Cycles and isotopes are explicitly unsupported',()=>{
 const model=PRESETS.ezbutE(); const bond=model.bonds.find(b=>b.order===2);
 model.atoms[0].isotope=13;assert.equal(analyzeEZ(model,bond.id).status,'unsupported');delete model.atoms[0].isotope;
 const outer=model.atoms.filter(a=>a.symbol==='C').map(a=>a.id);
 model.bonds.push({id:100,a:outer[0],b:outer[outer.length-1],order:1});
 assert.equal(analyzeEZ(model,bond.id).status,'unsupported');
});
test('Geometry near 90 degrees around C=C is reported as degenerate rather than E/Z',()=>{
 const m=PRESETS.ezbutE(),b=m.bonds.find(x=>x.order===2);
 // Rotate side (b.b) 90 degrees around the double bond axis, keeping all attached
 // ligands on the second center consistent to represent a twisted double bond.
 const a=m.atoms.find(x=>x.id===b.a),c=m.atoms.find(x=>x.id===b.b);
 const axis=[c.x-a.x,c.y-a.y,c.z-a.z],mag=Math.hypot(...axis),u=axis.map(v=>v/mag);
 const side=m.bonds.filter(x=>x.a===c.id||x.b===c.id).filter(x=>x.id!==b.id).map(x=>x.a===c.id?x.b:x.a);
 for(const id of side){const atom=m.atoms.find(x=>x.id===id),v=[atom.x-c.x,atom.y-c.y,atom.z-c.z];
  const cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]],parallel=u.reduce((s,t,i)=>s+t*v[i],0);
  [atom.x,atom.y,atom.z]=cross.map((w,i)=>w+parallel*u[i]+[c.x,c.y,c.z][i]);
 }
 assert.equal(analyzeEZ(m,b.id).status,'degenerate');
 assert.throws(()=>invertEZ(m,b.id),/non assegnabile/);
});
test('Viewing rotation does not affect E/Z because it depends on coordinates',()=>{
 const m=PRESETS.ezbutZ(),b=m.bonds.find(x=>x.order===2);
 const rotated=clone(m);
 for(const a of rotated.atoms){const x=a.x,y=a.y,z=a.z;
  a.x=x*Math.cos(.8)-y*Math.sin(.8);
  a.y=x*Math.sin(.8)+y*Math.cos(.8);
  a.z=z;
 }
 assert.equal(analyzeEZ(rotated,b.id).descriptor,'Z');
});
test('Mirror a molecule preserves E/Z (unlike R/S)',()=>{
 const m=PRESETS.ezbutE(),b=m.bonds.find(x=>x.order===2);
 mirrorMolecule(m);assert.equal(analyzeEZ(m,b.id).descriptor,'E');
});
test('Analysis is read-only and inversion is rejected atomically when unsupported',()=>{
 const m=PRESETS.ethene(),b=m.bonds.find(x=>x.order===2),snap=JSON.stringify(m);
 analyzeEZ(m,b.id);assert.equal(JSON.stringify(m),snap);
 assert.throws(()=>invertEZ(m,b.id));assert.equal(JSON.stringify(m),snap);
});

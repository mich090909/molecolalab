import test from 'node:test';
import assert from 'node:assert/strict';
import {PRESETS,blank,clone,findAtom,addAtom,setBond,dist,angle,dihedral,exportProject,importProject,exportCJSON,importCJSON} from '../src/core.js';
import {analyzeCIP,summarizeCIP,mirrorMolecule} from '../src/stereo.js';

/** Reference descriptors from RDKit 2025.09.4: ETKDG (seed 55), MMFF94, AssignStereochemistryFrom3D.
 * RDKit is used only for independent test case generation; it is not needed at runtime.
 */
const cases={brclfr:[[2,'R']],brclfs:[[2,'S']],butanolr:[[2,'R']],butanols:[[2,'S']],lactic_s:[[2,'S']],alanine_s:[[2,'S']],ditartaric:[[4,'R'],[6,'R']]};
const getAssigned=d=>summarizeCIP(d).filter(s=>s.status==='assigned').map(s=>[s.id,s.descriptor]);
const points=d=>d.atoms.map(a=>({...a}));
const lengths=d=>d.bonds.map(b=>dist(findAtom(d,b.a),findAtom(d,b.b)));
for(const [key,expected] of Object.entries(cases)){
 test(`CIP stereochemistry verified with independent RDKit reference: ${key}`,()=>{
  const d=PRESETS[key]();assert.deepEqual(getAssigned(d),expected);
  for(const x of summarizeCIP(d))assert.equal(x.status,'assigned');
 });
 test(`Mirroring reverses every R/S center without modifying bond lengths: ${key}`,()=>{
  const d=PRESETS[key](),before=lengths(d),original=getAssigned(d);
  mirrorMolecule(d);
  assert.deepEqual(getAssigned(d),original.map(([id,s])=>[id,s==='R'?'S':'R']));
  for(let k=0;k<before.length;k++)assert.ok(Math.abs(before[k]-lengths(d)[k])<1e-9);
  assert.ok(d.name.startsWith('Speculare di'));
  mirrorMolecule(d);assert.deepEqual(getAssigned(d),original);
 });
 test(`All CIP descriptors survive rigid 3D rotation and translation: ${key}`,()=>{
  const d=PRESETS[key](),original=getAssigned(d),theta=1.023,beta=-0.47;
  for(const a of d.atoms){
   const x=a.x,y=a.y,z=a.z;
   const x1=x*Math.cos(theta)-y*Math.sin(theta);
   const y1=x*Math.sin(theta)+y*Math.cos(theta);
   a.x=x1+12.4;a.y=y1*Math.cos(beta)-z*Math.sin(beta)-4.3;
   a.z=y1*Math.sin(beta)+z*Math.cos(beta)+0.25;
  }
  assert.deepEqual(getAssigned(d),original);
 });
 test(`R/S survives lossless CJSON and native JSON export/import: ${key}`,()=>{
  const d=PRESETS[key](),original=getAssigned(d);
  assert.deepEqual(getAssigned(importCJSON(exportCJSON(d))),original);
  assert.deepEqual(getAssigned(importProject(exportProject(d))),original);
 });
}
test('Highest atomic-number priority: Br > Cl > F > H, independent of coordinate order',()=>{
 const d=PRESETS.brclfr(),s=analyzeCIP(d,2);
 assert.deepEqual(s.priorities.map(x=>x.symbol),['Br','Cl','F','H']);
 assert.equal(s.descriptor,'R');
});
test('Secondary CIP rankings distinguish carboxyl from methyl in (S)-lactic acid',()=>{
 const d=PRESETS.lactic_s(),s=analyzeCIP(d,2);
 assert.deepEqual(s.priorities.map(x=>x.symbol),['O','C','C','H']);
 assert.equal(s.descriptor,'S');
 const neighbors=s.priorities.map(e=>e.atomId);
 assert.equal(neighbors[1],4); // carbonyl carbon precedes methyl
 assert.equal(neighbors[2],1);
});
test('CIP ranking of alanine: N > carboxyl C > methyl C > H',()=>{
 const d=PRESETS.alanine_s(),s=analyzeCIP(d,2);
 assert.deepEqual(s.priorities.map(x=>x.symbol),['N','C','C','H']);
 assert.equal(s.priorities[1].atomId,4);
 assert.equal(s.descriptor,'S');
});
test('Identical H atoms do not become fictitious stereocenters',()=>{
 for(const key of ['methane','ammonium','ethane','butane']){
  const d=PRESETS[key]();assert.equal(getAssigned(d).length,0,key);
  assert.equal(summarizeCIP(d).length,0,key);
 }
});
test('Unsuitable centers with implicit hydrogen are not assigned',()=>{
 const d=PRESETS.brclfr();d.bonds=d.bonds.filter(b=>b.a!==2||b.b!==5);d.atoms=d.atoms.filter(a=>a.id!==5);
 assert.equal(analyzeCIP(d,2).status,'not_tetrahedral');
});
test('Non tetrahedral nitrogen lone pair does not produce permanent R/S',()=>{
 const d=PRESETS.ammonia();assert.equal(analyzeCIP(d,1).status,'not_tetrahedral');
});
test('Cyclic side groups are conservatively unsupported instead of assigned',()=>{
 const d=PRESETS.brclfr();const a=addAtom(d,'C',4,1,0),b=addAtom(d,'C',5,2,0),c=addAtom(d,'C',3.5,2,0);
 setBond(d,3,a);setBond(d,a,b);setBond(d,b,c);setBond(d,c,a);
 assert.equal(analyzeCIP(d,2).status,'unsupported');
 assert.match(analyzeCIP(d,2).reason,/anello/i);
});
test('Isotopic substitution is not yet supported',()=>{
 const d=PRESETS.brclfr();findAtom(d,5).isotope=2;
 assert.equal(analyzeCIP(d,2).status,'unsupported');
});
test('Coplanar ligand coordinates do not generate a false R/S descriptor',()=>{
 const d=PRESETS.brclfr();for(const a of d.atoms)a.z=0;
 assert.equal(analyzeCIP(d,2).status,'degenerate');
});
test('Overlapping ligands do not yield R/S',()=>{
 const d=PRESETS.brclfr();findAtom(d,3).x=findAtom(d,4).x;findAtom(d,3).y=findAtom(d,4).y;findAtom(d,3).z=findAtom(d,4).z;
 assert.equal(analyzeCIP(d,2).status,'degenerate');
});
test('Swapping two distinct spatial ligands flips R/S but not the priority ranking',()=>{
 const d=PRESETS.brclfr(),a=findAtom(d,3),b=findAtom(d,4),xyz=[a.x,a.y,a.z];
 const p=analyzeCIP(d,2);[a.x,a.y,a.z]=[b.x,b.y,b.z];[b.x,b.y,b.z]=xyz;
 const q=analyzeCIP(d,2);assert.equal(p.descriptor,'R');assert.equal(q.descriptor,'S');
 assert.deepEqual(q.priorities,p.priorities);
});
test('Changing a ligand species triggers priority recalculation or unresolved chirality',()=>{
 const d=PRESETS.brclfr();findAtom(d,4).symbol='Cl';
 assert.notEqual(analyzeCIP(d,2).status,'assigned');
});
test('No molecule and unknown atoms handled safely',()=>{
 const d=blank();assert.deepEqual(summarizeCIP(d),[]);
 assert.equal(analyzeCIP(d,999).status,'unavailable');
 assert.throws(()=>mirrorMolecule(d),/atomi/i);
});
test('Mirror operation leaves net charges, all angles and all bond lengths unchanged',()=>{
 const d=PRESETS.lactic_s();
 const old=points(d);const initial=lengths(d);
 const oldAngle=angle(old[1],old[0],old[6]);
 mirrorMolecule(d);
 for(let i=0;i<initial.length;i++)assert.ok(Math.abs(initial[i]-lengths(d)[i])<1e-9);
 assert.equal(d.atoms.length,old.length);
 assert.equal(d.bonds.length,PRESETS.lactic_s().bonds.length);
 for(let i=0;i<old.length;i++)assert.equal(d.atoms[i].charge,old[i].charge);
});

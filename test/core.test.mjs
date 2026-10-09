import test from 'node:test';import assert from 'node:assert/strict';
import {setAtomPosition,setAtomElement,setBondLength,estimateVSEPR,netCharge,SYMBOLS,ELEMENTS,blank,clone,element,addAtom,setBond,eraseAtom,bondBetween,dist,angle,dihedral,formula,
 hydrogenate,autoBond,regularize,exportXYZ,importXYZ,exportMOL,importMOL,exportProject,importProject,PRESETS} from '../src/core.js';
const near=(v,expected,tol=.02)=>assert.ok(Math.abs(v-expected)<tol,`${v} non vicino a ${expected}`);
test('118 elementi, ognuno con Z progressivo',()=>{assert.equal(SYMBOLS.length,118);assert.equal(ELEMENTS.at(-1).symbol,'Og');assert.equal(element('C').Z,6);});
test('molecole precaricate coerenti per formula',()=>{for(const [key,expected] of Object.entries({water:'H2O',methane:'CH4',ammonia:'H3N',co2:'CO2',ethene:'C2H4',ethyne:'C2H2',formaldehyde:'CH2O',benzene:'C6H6',bf3:'BF3',sf6:'F6S'}))assert.equal(formula(PRESETS[key]()),expected,key);});
test('geometria H2O: lunghezze OH e angolo H-O-H',()=>{const d=PRESETS.water();near(dist(d.atoms[0],d.atoms[1]),.957,.005);near(angle(d.atoms[1],d.atoms[0],d.atoms[2]),104.5,.2);});
test('metano tetraedrico',()=>{const d=PRESETS.methane();near(angle(d.atoms[1],d.atoms[0],d.atoms[2]),109.471,.01);});
test('CO2 lineare',()=>{const d=PRESETS.co2();near(angle(d.atoms[0],d.atoms[1],d.atoms[2]),180,.001);});
test('distanze, angoli, diedri: casi semplici',()=>{const a={x:1,y:0,z:0},b={x:0,y:0,z:0},c={x:0,y:1,z:0},d={x:0,y:1,z:1};near(dist(a,b),1);near(angle(a,b,c),90);near(Math.abs(dihedral(a,b,c,d)),90);});
test('aggiunta e modifica ordine legame, unicità',()=>{const d=blank();const a=addAtom(d,'C',0,0,0),b=addAtom(d,'O',1.2,0,0);setBond(d,a,b,2);assert.equal(d.bonds.length,1);setBond(d,b,a,3);assert.equal(d.bonds.length,1);assert.equal(bondBetween(d,a,b).order,3);});
test('atomo sconosciuto respinto',()=>{assert.throws(()=>addAtom(blank(),'Xx',0,0,0));});
test('elimina atomo e legami',()=>{const d=PRESETS.water();eraseAtom(d,1);assert.equal(d.atoms.length,2);assert.equal(d.bonds.length,0);});
test('completa H su carbonio isolato: CH4',()=>{const d=blank();addAtom(d,'C',0,0,0);assert.equal(hydrogenate(d),4);assert.equal(formula(d),'CH4');assert.equal(hydrogenate(d),0);});
test('completamento valenza semplice del metanale',()=>{const d=blank();addAtom(d,'C',0,0,0);addAtom(d,'O',1.2,0,0,1,2);const n=hydrogenate(d);assert.equal(n,2);assert.equal(formula(d),'CH2O');});
test('non aggiunge H a BF3 o SF6',()=>{for(const name of ['bf3','sf6'])assert.equal(hydrogenate(PRESETS[name]()),0,name);});
test('riconosce legami XYZ con distanza vicina e nessuno a distanza grande',()=>{const d=blank();addAtom(d,'H',0,0,0);addAtom(d,'H',.75,0,0);addAtom(d,'He',8,0,0);assert.equal(autoBond(d),1);assert.equal(d.bonds[0].order,1);});
test('XYZ round-trip coordinate e formule',()=>{for(const name of ['water','methane','benzene']){const d=PRESETS[name](),r=importXYZ(exportXYZ(d));assert.equal(r.atoms.length,d.atoms.length);assert.equal(formula(r),formula(d));near(r.atoms[0].x,d.atoms[0].x,.00001);}});
test('MOL V2000 round-trip con ordini di legame',()=>{for(const name of ['ethene','ethyne','benzene','bf3']){const d=PRESETS[name](),r=importMOL(exportMOL(d));assert.equal(r.bonds.length,d.bonds.length);assert.deepEqual(r.bonds.map(b=>b.order),d.bonds.map(b=>b.order));assert.equal(formula(r),formula(d));}});
test('MOL conserva carica formale',()=>{const d=blank();addAtom(d,'N',0,0,0);d.atoms[0].charge=1;const r=importMOL(exportMOL(d));assert.equal(r.atoms[0].charge,1);});
test('JSON nativo round-trip identico',()=>{const d=PRESETS.methane();assert.deepEqual(importProject(exportProject(d)),d);});
test('protezione da formato JSON non proprio',()=>{assert.throws(()=>importProject('{"unrelated":123}'));});
test('protegge da coordinate non numeriche',()=>{assert.throws(()=>importXYZ('1\nTest\nC nan 0 0'));});
test('riordino geometrico converge senza NaN',()=>{const d=PRESETS.water();regularize(d,100);for(const a of d.atoms)assert.ok([a.x,a.y,a.z].every(Number.isFinite));});
test('clone non condivide atomi con originale',()=>{const d=PRESETS.water();const c=clone(d);c.atoms[0].x=100;assert.notEqual(d.atoms[0].x,c.atoms[0].x);});

// Versione 0.2 — geometric editing, charge, VSEPR diagnostics and import integrity.
test('nuovi esempi: formula e cariche coerenti',()=>{
 for(const [key,f,q] of [['pcl5','Cl5P',0],['hcn','CHN',0],['ammonium','H4N',1],['carbonate','CO3',-2]]){
  const d=PRESETS[key]();assert.equal(formula(d),f,key);assert.equal(netCharge(d),q,key);
 }
});
test('VSEPR: acqua angolare AX2E2, ammoniaca piramidale AX3E',()=>{
 const water=estimateVSEPR(PRESETS.water(),1),ammonia=estimateVSEPR(PRESETS.ammonia(),1);
 assert.equal(water.formula,'AX2E2');assert.equal(water.molecularGeometry,'angolare');
 assert.equal(water.electronicGeometry,'tetraedrica');
 assert.equal(ammonia.formula,'AX3E');assert.equal(ammonia.molecularGeometry,'piramidale trigonale');
});
test('VSEPR: metano tetraedrico, CO2 lineare, BF3 trigonale',()=>{
 for(const [key,center,shape,AXE] of [['methane',1,'tetraedrica','AX4'],['co2',2,'lineare','AX2'],['bf3',1,'trigonale planare','AX3']]){
  const v=estimateVSEPR(PRESETS[key](),center);assert.ok(v.supported,key);assert.equal(v.formula,AXE);assert.equal(v.molecularGeometry,shape);
 }
});
test('VSEPR: SF6, PCl5 e ammonio (catione) corretti',()=>{
 for(const [key,shape,AXE] of [['sf6','ottaedrica','AX6'],['pcl5','bipiramidale trigonale','AX5'],['ammonium','tetraedrica','AX4']]){
  const v=estimateVSEPR(PRESETS[key](),1);assert.ok(v.supported,key);assert.equal(v.molecularGeometry,shape);assert.equal(v.formula,AXE);
 }
});
test('VSEPR: carbonato, singoli atomi e metalli non presumono geometrie non fondate',()=>{
 const d=PRESETS.carbonate();assert.equal(estimateVSEPR(d,1).molecularGeometry,'trigonale planare');
 assert.equal(estimateVSEPR(d,2).supported,false);
 const x=blank();addAtom(x,'Fe',0,0,0);addAtom(x,'H',1,0,0,1);addAtom(x,'H',-1,0,0,1);
 assert.equal(estimateVSEPR(x,1).supported,false);
});
test('completamento H non neutralizza anioni carbonato',()=>{
 const d=PRESETS.carbonate();assert.equal(hydrogenate(d),0);assert.equal(formula(d),'CO3');assert.equal(netCharge(d),-2);
});
test('modifica elemento e coordinate valida con protezione limiti',()=>{
 const d=PRESETS.methane();setAtomElement(d,1,'N');d.atoms[0].charge=1;
 assert.equal(formula(d),'H4N');assert.equal(netCharge(d),1);
 setAtomPosition(d,1,.5,-.2,1.25);assert.deepEqual([d.atoms[0].x,d.atoms[0].y,d.atoms[0].z],[.5,-.2,1.25]);
 assert.throws(()=>setAtomPosition(d,1,Infinity,0,0));assert.throws(()=>setAtomPosition(d,1,1001,0,0));
 assert.throws(()=>setAtomElement(d,1,'Xx'));assert.equal(d.atoms[0].symbol,'N');
});
test('modifica lunghezza: sposta solo il secondo estremo e mantiene legami',()=>{
 const d=PRESETS.co2(),bond=d.bonds[0],fixed={...d.atoms[0]};
 const centerBefore={...d.atoms[1]};setBondLength(d,bond.id,1.50);
 near(dist(d.atoms[0],d.atoms[1]),1.5,.00001);
 assert.deepEqual(d.atoms[0],fixed);assert.deepEqual(d.atoms[2],PRESETS.co2().atoms[2]);
 assert.notDeepEqual(d.atoms[1],centerBefore);
 assert.throws(()=>setBondLength(d,bond.id,.15));
 assert.throws(()=>setBondLength(d,999,1));
});
test('inserimento atomo fallito lascia il documento intatto',()=>{
 const d=blank();addAtom(d,'C',0,0,0);const before=clone(d);
 assert.throws(()=>addAtom(d,'O',1,0,0,99,1));assert.deepEqual(d,before);
 assert.throws(()=>addAtom(d,'O',1,0,0,1,4));assert.deepEqual(d,before);
});
test('import JSON rifiuta id legami duplicati, legami ripetuti e cariche invalide',()=>{
 const d=PRESETS.water();const x=clone(d);x.bonds[1].id=x.bonds[0].id;
 assert.throws(()=>importProject(JSON.stringify({format:'molecolalab-v1',molecule:x})));
 const y=clone(d);y.bonds[1].a=1;y.bonds[1].b=2;
 assert.throws(()=>importProject(JSON.stringify({format:'molecolalab-v1',molecule:y})));
 const z=clone(d);z.atoms[0].charge=99;
 assert.throws(()=>importProject(JSON.stringify({format:'molecolalab-v1',molecule:z})));
});
test('nuove molecole mantengono legami e cariche nei formati supportati',()=>{
 for(const key of ['pcl5','hcn','ammonium','carbonate']){
  const d=PRESETS[key](),r=importMOL(exportMOL(d));
  assert.equal(formula(r),formula(d));assert.equal(netCharge(r),netCharge(d));
  assert.deepEqual(r.bonds.map(b=>b.order),d.bonds.map(b=>b.order));
 }
});

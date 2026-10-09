import test from 'node:test';
import assert from 'node:assert/strict';
import {
 PRESETS,blank,clone,addAtom,setBond,findAtom,dist,dihedral,formula,netCharge,
 setDihedral,exportCJSON,importCJSON,unhandledCJSONFields
} from '../src/core.js';
const almost=(actual,expected,epsilon=1e-4)=>assert.ok(Math.abs(actual-expected)<epsilon,`ottenuto ${actual}, previsto ${expected}`);
const coords=d=>d.atoms.map(a=>[a.x,a.y,a.z]);
const torsion=d=>dihedral(...[1,2,3,4].map(id=>findAtom(d,id)));
const pairs=d=>d.bonds.map(b=>dist(findAtom(d,b.a),findAtom(d,b.b)));
const allPresets=Object.keys(PRESETS);

test('CJSON: esportazione secondo la struttura Chemical JSON v1 e indici zero-based',()=>{
 const d=PRESETS.water(),x=JSON.parse(exportCJSON(d));
 assert.equal(x.chemicalJson,1);assert.equal(x.atoms.elements.number[0],8);
 assert.deepEqual(x.atoms.elements.number,[8,1,1]);
 assert.equal(x.atoms.coords['3d'].length,9);
 assert.deepEqual(x.bonds.connections.index,[0,1,0,2]);
 assert.deepEqual(x.bonds.order,[1,1]);
 assert.deepEqual(x.atoms.formalCharges,[0,0,0]);
});
test('CJSON: round-trip atomi, ordine, coordinate, cariche per tutti i preset',()=>{
 for(const key of allPresets){
  const d=PRESETS[key](),r=importCJSON(exportCJSON(d));
  assert.deepEqual(r,d,`roundtrip ${key}`);
 }
});
test('CJSON: preserva correttamente legami multipli e cariche anioniche',()=>{
 const carbonate=JSON.parse(exportCJSON(PRESETS.carbonate()));
 assert.deepEqual(carbonate.bonds.order,[2,1,1]);
 assert.equal(carbonate.atoms.formalCharges.filter(x=>x===-1).length,2);
 const back=importCJSON(JSON.stringify(carbonate));assert.equal(netCharge(back),-2);
 assert.equal(formula(back),'CO3');
});
test('CJSON: import senza bonds non inventa legami',()=>{
 const x=JSON.parse(exportCJSON(PRESETS.water()));delete x.bonds;
 const r=importCJSON(JSON.stringify(x));assert.equal(r.atoms.length,3);assert.equal(r.bonds.length,0);
});
test('CJSON: import senza formalCharges assume atomi neutri',()=>{
 const x=JSON.parse(exportCJSON(PRESETS.water()));delete x.atoms.formalCharges;
 const r=importCJSON(JSON.stringify(x));assert.equal(netCharge(r),0);
});
test('CJSON: respinge coordinate incomplete, non finite e fuori scala',()=>{
 const x=JSON.parse(exportCJSON(PRESETS.water()));x.atoms.coords['3d'].pop();
 assert.throws(()=>importCJSON(JSON.stringify(x)),/coordinate/i);
 x.atoms.coords['3d'].push(0);x.atoms.coords['3d'][0]=1001;
 assert.throws(()=>importCJSON(JSON.stringify(x)),/coordinate/i);
});
test('CJSON: respinge numeri atomici errati e cariche formalmente invalide',()=>{
 const x=JSON.parse(exportCJSON(PRESETS.water()));x.atoms.elements.number[1]=119;
 assert.throws(()=>importCJSON(JSON.stringify(x)),/numero atomico/i);
 x.atoms.elements.number[1]=1;x.atoms.formalCharges=[0,0];
 assert.throws(()=>importCJSON(JSON.stringify(x)),/cariche/i);
});
test('CJSON: respinge ordine aromatico 4, frazionario 1.5 e legami duplicati',()=>{
 const original=JSON.parse(exportCJSON(PRESETS.water()));
 for(const order of [4,1.5,0]){const x=clone(original);x.bonds.order[0]=order;assert.throws(()=>importCJSON(JSON.stringify(x)),/legame/i);}
 const copy=clone(original);copy.bonds.order.push(1);copy.bonds.connections.index.push(1,0);
 assert.throws(()=>importCJSON(JSON.stringify(copy)),/duplicato/i);
});
test('CJSON: respinge indici fuori intervallo e connettività non coerente',()=>{
 const original=JSON.parse(exportCJSON(PRESETS.water()));
 for(const [a,b] of [[1,1],[0,6],[-1,2]]){
  const x=clone(original);x.bonds.connections.index.splice(0,2,a,b);
  assert.throws(()=>importCJSON(JSON.stringify(x)),/legame/i);
 }
});
test('CJSON: previene perdita silenziosa di dati di carica totale non trasferibile',()=>{
 const x=JSON.parse(exportCJSON(PRESETS.water()));x.properties={totalCharge:1};
 assert.throws(()=>importCJSON(JSON.stringify(x)),/carica totale/i);
});
test('CJSON: rifiuta celle periodiche non gestite',()=>{
 const x=JSON.parse(exportCJSON(PRESETS.water()));x.unitCell={a:2,b:2,c:2,alpha:90,beta:90,gamma:90};
 assert.throws(()=>importCJSON(JSON.stringify(x)),/periodiche/i);
});
test('CJSON: rileva dati aggiuntivi, ma non segnala il totale carica supportato',()=>{
 const x=JSON.parse(exportCJSON(PRESETS.ammonium()));x.properties={totalCharge:1};
 assert.deepEqual(unhandledCJSONFields(x),[]);
 x.vibrations={frequencies:[1000]};x.atoms.labels=['N','H','H','H','H'];
 assert.deepEqual(unhandledCJSONFields(x).sort(),['atoms.labels','vibrations']);
});
test('CJSON: documenti diversi rifiutati senza sostituzione del precedente',()=>{
 assert.throws(()=>importCJSON('non JSON'),/JSON non valido/i);
 assert.throws(()=>importCJSON('{}'),/chemicalJson/i);
 assert.throws(()=>importCJSON('{"chemicalJson":0,"atoms":{}}'),/chemicalJson/i);
});

test('torsione: atomi del butano e legami MMFF94 di riferimento',()=>{
 const d=PRESETS.butane();assert.equal(d.atoms.length,14);assert.equal(d.bonds.length,13);
 for(const length of pairs(d).slice(0,3))assert.ok(length>1.45&&length<1.57,length);
 assert.ok(Number.isFinite(torsion(d)));
});
test('torsione: imposta 60°, −60°, 180°, −180°, 0° senza accumulare errore',()=>{
 const d=PRESETS.butane();
 for(const value of [60,-60,180,-180,0,92.5]){
  const res=setDihedral(d,1,2,3,4,value);assert.ok(res.moved>=1);
  const achieved=torsion(d);almost(Math.abs(((achieved-value+540)%360)-180),0,1e-6);
 }
});
test('torsione: resta intatta la parte fissa e gli atomi sull’asse',()=>{
 const d=PRESETS.butane(),before=clone(d);setDihedral(d,1,2,3,4,175);
 for(const id of [1,2,3,5,6,7,8,9])assert.deepEqual(findAtom(d,id),findAtom(before,id));
 assert.notDeepEqual(findAtom(d,4),findAtom(before,4));
});
test('torsione: tutte le distanze dei legami sono invarianti',()=>{
 const d=PRESETS.butane(),initial=pairs(d);
 setDihedral(d,1,2,3,4,-160);
 for(let i=0;i<initial.length;i++)almost(pairs(d)[i],initial[i],1e-9);
});
test('torsione: non muove altri frammenti scollegati',()=>{
 const d=PRESETS.butane();const id=addAtom(d,'He',20,0,0);const original={...findAtom(d,id)};
 setDihedral(d,1,2,3,4,120);assert.deepEqual(findAtom(d,id),original);
});
test('torsione: doppio legame centrale bloccato',()=>{
 const d=PRESETS.butane();setBond(d,2,3,2);const before=clone(d);
 assert.throws(()=>setDihedral(d,1,2,3,4,30),/legami centrali singoli/i);
 assert.deepEqual(d,before);
});
test('torsione: anello benzene bloccato anche se un legame è fatto singolo',()=>{
 const d=PRESETS.benzene();const before=clone(d);
 // 1–6 is one of the alternating single bonds (6-1 originally 1). Ring stays closed.
 assert.throws(()=>setDihedral(d,2,1,6,5,30),/anello/i);
 assert.deepEqual(d,before);
});
test('torsione: legami non consecutivi, ripetizioni, atomi inesistenti bloccati',()=>{
 const d=PRESETS.butane();const before=clone(d);
 assert.throws(()=>setDihedral(d,1,2,4,3,30),/legame centrale/i);
 assert.throws(()=>setDihedral(d,1,2,2,3,30),/distinti/i);
 assert.throws(()=>setDihedral(d,1,2,3,99,30),/non esistono/i);
 assert.deepEqual(d,before);
});
test('torsione: angoli fuori [−180,180] o NaN bloccati',()=>{
 const d=PRESETS.butane();const before=clone(d);
 for(const x of [181,-181,NaN,Infinity])assert.throws(()=>setDihedral(d,1,2,3,4,x),/richiesto/i);
 assert.deepEqual(d,before);
});
test('torsione: coordinate allineate e diedro indefinito bloccati',()=>{
 const d=blank();let a=addAtom(d,'C',0,0,0),b=addAtom(d,'C',1,0,0),c=addAtom(d,'C',2,0,0),d4=addAtom(d,'C',2,1,0);
 setBond(d,a,b);setBond(d,b,c);setBond(d,c,d4);
 assert.throws(()=>setDihedral(d,a,b,c,d4,90),/Diedro non definito/i);
});

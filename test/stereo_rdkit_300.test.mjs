import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {addAtom,blank,findAtom,setBond} from '../src/core.js';
import {analyzeCIP} from '../src/stereo.js';

/** Fixed, independently labeled reference structures generated offline with
 * RDKit 2025.09.4, ETKDG 3D embed (seed 55), MMFF94 and stereochemistry from 3D.
 * All snapshots are standard acyclic organic structures and have explicit H.
 * The JavaScript app has zero dependency on RDKit.
 */
const cases=JSON.parse(readFileSync(new URL('./fixtures/cip_rdkit_300.json',import.meta.url),'utf8'));
assert.equal(cases.length,300,'Expected exactly 300 independent controls');
for(const [i,c] of cases.entries()){
 test(`CIP reference RDKit ${String(i+1).padStart(3,'0')} — ${c.smi}`,()=>{
  const doc=blank();
  for(const [symbol,x,y,z,charge] of c.atoms){
   const id=addAtom(doc,symbol,x,y,z);
   findAtom(doc,id).charge=charge;
  }
  for(const [a,b,order] of c.bonds)setBond(doc,a,b,order);
  for(const [id,reference] of c.expected){
   const outcome=analyzeCIP(doc,id);
   assert.equal(outcome.status,'assigned',`${c.smi}: ${outcome.reason??'no result'}`);
   assert.equal(outcome.descriptor,reference,c.smi);
  }
 });
}

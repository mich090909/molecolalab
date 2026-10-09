/** Rebuild the offline demo from the actual codebase, with no external dependencies. */
import {readFileSync,writeFileSync} from 'node:fs';
import {resolve} from 'node:path';
const root=resolve(import.meta.dirname,'..');
const read=name=>readFileSync(resolve(root,name),'utf8');
const html=read('index.html'),css=read('src/style.css');
const core=read('src/core.js').replace(/^export /gm,'');
const stereo=read('src/stereo.js').replace(/^import\s*\{[\s\S]*?\}\s*from\s*['"]\.\/core\.js['"];?\s*/m,'').replace(/^export /gm,'');
const app=read('src/app.js');
const pattern=/^import\s*\{[\s\S]*?\}\s*from\s*['"]\.\/core\.js['"];?\s*/;
if(!pattern.test(app))throw Error('Struttura import app.js inattesa');
const stereoImport=/^import\s*\{[\s\S]*?\}\s*from\s*['"]\.\/stereo\.js['"];?\s*/m;
if(!stereoImport.test(app))throw Error('Import di stereochimica mancante');
const inlineApp=app.replace(pattern,'').replace(stereoImport,'');
if(!html.includes('<link rel="stylesheet" href="./src/style.css" />')||!html.includes('<script type="module" src="./src/app.js"></script>'))throw Error('Punti di sostituzione non trovati');
const combined=html.replace('<link rel="stylesheet" href="./src/style.css" />',`<style>\n${css}\n</style>`)
 .replace('<script type="module" src="./src/app.js"></script>',`<script type="module">\n${core}\n${stereo}\n${inlineApp}\n</script>`);
writeFileSync(resolve(root,'standalone.html'),combined);
process.stdout.write(`standalone.html generato (${combined.length} caratteri)\n`);

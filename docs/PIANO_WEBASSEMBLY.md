# Piano di integrazione del vero motore Avogadro WebAssembly

**Verifica effettuata il 9 ottobre 2026.** Il [PR OpenChemistry/avogadroapp#876](https://github.com/OpenChemistry/avogadroapp/pull/876) è *open/draft* e il repository [OpenChemistry/avogadro-web](https://github.com/OpenChemistry/avogadro-web) non distribuisce ancora un pacchetto installabile chiavi in mano. Il PR propone un editor Qt WebAssembly con controlli HTML, plugin Avogadro, parser/scrittori e un controllo JavaScript.

## Vincoli tecnici

- Qt 6.11.1 WebAssembly **multithread** ed Emscripten 4.0.7 secondo il README sperimentale della proposta.
- Compilazione di dipendenze AvogadroLibs con opzioni `AVOGADRO_WEB_EDITOR=ON`, sottoprocessi e Python disabilitati; installazione in directory pulite.
- Non basta la pagina HTML: occorrono `avogadro-web.wasm`, loader Qt, altri asset e hosting con `Cross-Origin-Opener-Policy: same-origin`, `Cross-Origin-Embedder-Policy: require-corp` (ed eventualmente `Cross-Origin-Resource-Policy: same-origin`). La variante multithread richiede anche WebGL2, SharedArrayBuffer e HTTPS o localhost.
- **GitHub Pages** non offre normalmente header personalizzati per la cross-origin isolation. Hosting alternativo con header controllabili o diversa build possono essere necessari.
- Verifica delle licenze Qt LGPLv3/GPLv3, delle limitazioni del linking statico WebAssembly e delle licenze di ogni dipendenza/plugin, prima di ridistribuire file compilati; rispettare le condizioni di relinking, source e notices.
- La fase iniziale non include notebook, SMILES, plugin Python, calcoli server, tablet gestures e chimica quantistica completa: non promettere parità funzionale con Avogadro desktop.

## Percorso incrementale e reversibile

1. **Stabile v0.2:** mantenere un riferimento Git immutabile (branch/tag) prima di aggiornare il prodotto principale.
2. **Interoperabilità v0.3:** lettura e scrittura CJSON per elementi, coordinate, cariche e legami. Funzioni geometriche locali come torsioni modificabili. Nessun binario Avogadro presente.
3. **Sperimentazione separata:** verificare licenze e creare *build riproducibile* del PR in ambiente SDK compatibile, solo in branch isolato. Non copiare un binario da una demo non verificata.
4. **Adapter:** testare trasferimento in entrambe le direzioni tramite CJSON/MOL e confrontare topologia, cariche, coordinate e un set di molecole note.
5. **Validazione:** test su desktop, Chrome e Safari/iPad, con fallback al motore indipendente esistente; raccogliere prestazioni, memory footprint e licenze nel pacchetto.
6. **Eventuale pubblicazione:** decisione separata dell'utente, con permessi e hosting esplicitamente approvati.

Fonti principali: [PR WebAssembly](https://github.com/OpenChemistry/avogadroapp/pull/876), [Qt for WebAssembly](https://doc.qt.io/qt-6/wasm.html), [Qt LGPL obligations](https://www.qt.io/development/open-source-lgpl-obligations), [Chemical JSON](https://github.com/OpenChemistry/chemicaljson).

# Licenze, copyright e limiti scientifici

**Audit preliminare, 9 ottobre 2026. Non costituisce parere legale.**

## Stato di distribuzione, v0.1

Questa versione usa esclusivamente codice JavaScript/HTML/CSS originale scritto per MolecolaLab Web. Non include binari Qt, codice copiato da Avogadro, il logo Avogadro né librerie esterne. Licenza applicata al codice originale: MIT, con testo in `/LICENSE`.

Riferimenti a software terzo non significano inclusione dei relativi sorgenti o binari.

## Progetti terzi per la futura integrazione

| Dipendenza | Stato | Licenza/attenzione | Riferimento |
| --- | --- | --- | --- |
| Avogadro 2 `avogadroapp` | Non inclusa | BSD 3 clausole; conservare attribuzioni, testo e disclaimer | https://github.com/OpenChemistry/avogadroapp/blob/master/LICENSE |
| Avogadro `avogadrolibs` | Non inclusa | BSD 3 clausole per il progetto; verificare moduli e tutte le dipendenze | https://github.com/OpenChemistry/avogadrolibs/blob/master/LICENSE |
| `avogadro-web` | Non inclusa | BSD 3 clausole; progetto iniziale | https://github.com/OpenChemistry/avogadro-web |
| PR `avogadroapp#876` | Non inclusa | Implementazione prototipale non ancora integrata al 9 ottobre 2026; verificare revisioni | https://github.com/OpenChemistry/avogadroapp/pull/876 |
| Qt per WebAssembly | Non inclusa | LGPLv3 e/o GPLv3 a seconda dei moduli/versioni; particolare cautela con linking statico e obblighi di ri-link / ricompilazione | https://www.qt.io/development/open-source-lgpl-obligations |
| Plugin QTAIM/GPL | Non incluso | Non assumere che abbia licenza BSD | https://avogadro.cc/develop/build.html |

### Regole prima di una release con motore Avogadro

1. Fissare esattamente commit e versioni delle librerie usate.
2. Eseguire un inventario SPDX di **tutte** le dipendenze effettivamente compilate, plugin inclusi.
3. Conservare i file COPYRIGHT/NOTICE/LICENSE e includere la documentazione prevista per distribuzione binaria.
4. Verificare l'effettiva licenza dei moduli Qt e come adempiere agli obblighi di disponibilità dei sorgenti, modifica/ri-linking e istruzioni per ricompilare e usare versioni modificate in WASM.
5. Se la catena Qt WASM richiede GPLv3, rispettare le condizioni GPLv3 per il lavoro combinato o scegliere un'architettura diversa. Non distribuire finché il controllo non è completato.
6. Non usare loghi, nomi di dominio o messaggi che suggeriscano una sponsorizzazione degli autori originali senza autorizzazione.
7. Rivedere licenze e condizioni **a ogni aggiornamento** di versione o dipendenza.

## Limiti scientifici v0.1

- Il rendering si basa su proiezioni di coordinate cartesiane 3D in Canvas; non è un motore OpenGL/WebGL Avogadro.
- Il calcolo delle distanze (Å), degli angoli e diedri è geometrico e dipende dalle coordinate presenti, senza stime di incertezza.
- Il comando «Completa H» usa valenze *neutral-form* semplici e non modella adeguatamente specie cariche, radicali, composti di coordinazione, aromaticità o valenze non ordinarie.
- La funzione «Riordina spaziatura» (fino a 220 atomi) applica vincoli e repulsioni geometriche ad hoc: **non ottimizza l'energia** e non produce geometrie conformazionali validate.
- Import XYZ: solo coordinate e tipo atomico, con inferenza dei legami singoli da distanze. XYZ non rappresenta bond-order; non recupera automaticamente carica/risonanza.
- Import MOL V2000: usa legami di ordine 1–3; non gestisce completamente stereo, aromatizzazione, proprietà avanzate o SDF multi-record. L'importazione di SDF legge il primo record.
- Le viste «sfere di van der Waals» usano dimensioni illustrative e non vanno interpretate come dati quantitativi di raggi VdW.
- Il prototipo non contiene routine di calcolo di orbitali, meccanica molecolare, spettroscopia o dinamica.

## Privacy e dati

I file molecolari vengono letti localmente dal browser e non inviati a servizi remoti. Il sito non include analytics né telemetria. Attenzione: un futuro hosting di terze parti può conservare log HTTP ordinari, indipendenti dal codice dell'app.

## Aggiornamento v0.2 — 9 ottobre 2026

- Nessuna dipendenza esterna aggiunta: VSEPR semplificato, drag nel piano della vista, letture geometriche ed esportazioni restano codice originale MIT. **Non** è incorporato alcun binario/sorgente Avogadro, Qt o altro software terzo.
- Le nuove geometrie PCl₅ e carbonato sono idealizzazioni didattiche. PCl₅ è mostrato nella sua geometria molecolare bipiramidale trigonale, tipica della specie molecolare in fase gassosa; allo stato solido esistono altre specie ioniche.
- Il carbonato mostra una sola formula limite; risonanza e delocalizzazione non sono calcolate.
- VSEPR fornisce una **stima** dei domini su modelli Lewis a guscio chiuso e solo per alcuni elementi del blocco p; non è una validazione generale delle strutture, né gestisce in modo completo specie con legami complessi, radicali, transizioni elettroniche, metalli o ipervalenti inattesi.
- Trascinamento e cambio di lunghezza spostano coordinate, non effettuano minimizzazione energetica. Le misure sono esatte solo rispetto alle coordinate numeriche del modello; non sono risultati di ricerca.
- La saturazione H evita gli atomi con cariche formali diverse da zero: scelta conservativa per non neutralizzare, ad esempio, uno ione carbonato; altre valenze non ordinarie restano possibili e non validate.

## Aggiornamento v0.3: Chemical JSON e coordinate di riferimento

La versione 0.3 implementa **da zero** un convertitore del formato documentato [Chemical JSON](https://github.com/OpenChemistry/chemicaljson) secondo lo schema pubblico: nessun codice o binario del progetto Avogadro è stato incluso. Il formato può trasportare anche orbitali, spettri, dati periodicità, vibrazioni e proprietà varie, che non sono implementate nel modello di MolecolaLab. L'interfaccia rifiuta valori incompatibili (ad esempio legami aromatici non rappresentabili e celle periodiche) e richiede la conferma prima di scartare altre proprietà non supportate. Conservare sempre il sorgente CJSON completo.

Le coordinate 3D di **etano e butano** preimpostati sono state generate offline con RDKit 2025.09.4 mediante embedding tridimensionale e minimizzazione MMFF94. **RDKit non è distribuito**, né utilizzato dal browser, e la rotazione del diedro non è un'ottimizzazione MMFF/UFF. Documentazione RDKit: https://www.rdkit.org/ ; licenza repository https://github.com/rdkit/rdkit/blob/master/license.txt.

Le caratteristiche tecniche e gli obblighi di Qt sono documentati in `PIANO_WEBASSEMBLY.md`. Non pubblicare binari Avogadro WebAssembly senza aver verificato le condizioni di LGPL/GPL e ogni altra dipendenza.

## Aggiornamento v0.4 — stereochimica CIP didattica

Il modulo originale `src/stereo.js` è scritto specificamente per MolecolaLab, senza librerie Avogadro, Qt o RDKit incluse. È rilasciato con la medesima licenza MIT del progetto. Gli esempi 3D di molecole stereogeniche e **300 fixture** di confronto sono stati generati offline e confrontati con RDKit 2025.09.4 (embedding ETKDG, coordinate MMFF94, etichette 3D). I dati delle fixture non includono codice RDKit e l'applicazione non utilizza RDKit in esecuzione.

La v0.4 implementa solo un **sottoinsieme** delle regole CIP: centri tetraedrici C/Si/N⁺ con quattro sostituenti espliciti, strutture acicliche, ordinamenti per numero atomico e confronti ricorsivi (duplicazioni elementari su legami multipli). Non tratta anelli, isotopi, stereogenesi diversa dal centro tetraedrico, pseudoasimmetria, H impliciti, descrittori dipendenti da altre unità stereogeniche né regole CIP complete; casi fuori supporto non sono assegnati. Le coordinate rappresentano geometrie di esempio, non un risultato sperimentale. Documentazione: `STEREOCHIMICA_CIP.md`.

Alla versione 0.5, questo codice non incorpora binari WebAssembly di Avogadro né componenti Qt o RDKit. Per la distribuzione pubblica statica vale la licenza MIT indicata in `/LICENSE`. L’accessibilità e i requisiti normativi applicabili all’uso scolastico devono essere valutati separatamente; non è dichiarata conformità WCAG.


## Versione 0.5: E/Z e guida integrata

Il modulo E/Z è codice JavaScript originale distribuito sotto la licenza MIT del progetto. I modelli di riferimento sono stati generati offline mediante RDKit 2025.09.4, usato soltanto nello sviluppo per produrre coordinate e verificare 180 assegnazioni E/Z. **Il programma distribuito non include RDKit, Qt o Avogadro né incorpora i relativi componenti.** La generazione dell'altro isomero geometrico non è una torsione fisicamente libera del doppio legame. La guida completa è incorporata nel programma `standalone.html` per essere consultabile senza servizi esterni. Il confronto CIP resta limitato ai casi dichiarati in `docs/ISOMERIA_EZ.md`.

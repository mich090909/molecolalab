# MolecolaLab Web v0.5 — Modellazione molecolare e stereochimica R/S ed E/Z

Editor molecolare tridimensionale **didattico e indipendente**, gratuito, eseguibile nel browser **anche offline**. Realizzato in JavaScript/HTML/CSS originale e distribuito con licenza MIT. **Non include** codice o binari di Avogadro, Qt o RDKit.

## Stereochimica R/S (dalla versione 0.4)

- Individuazione e assegnazione **R/S** di centri tetraedrici supportati a partire da **connettività e coordinate 3D** (non dall'orientamento della telecamera).
- Priorità CIP da 1 a 4 con confronto ricorsivo dei sostituenti, inclusa la rappresentazione dei legami multipli di tipo carbonile e simili in strutture **acicliche**.
- Etichette R/S sui centri del modellino, numeri di priorità sugli atomi adiacenti quando selezioni il centro, pannello di spiegazione.
- Pulsante **Crea la molecola speculare**: riflette **tutti** i centri stereogenici; operazione reversibile con Annulla/Ripeti.
- Esempi: bromoclorofluorometano (R/S), butan-2-olo (R/S), (S)-acido lattico, (S)-alanina, (2R,3R)-acido tartarico.
- Geometrie degli esempi generate offline con RDKit 2025.09.4 (embedding ETKDG, MMFF94); **RDKit non è necessario** all'esecuzione né viene distribuito.

**Campo di applicazione:** la v0.5 assegna R/S a centri tetraedrici **C, Si e N⁺** con **quattro sostituenti rappresentati esplicitamente e quattro legami singoli**, in grafi molecolari aciclici. Non gestisce il confronto CIP completo IUPAC: **anelli, isotopi, idrogeni impliciti, pseudoasimmetria, descrittori dipendenti da altre unità stereogeniche e altre regole avanzate** non vengono assegnati. La connettività importata potrebbe anche non riflettere correttamente valenze e stereochimica. In tali casi il programma segnala l'impossibilità di concludere anziché attribuire R/S. **Non usare questa versione come strumento di certificazione scientifica.**

La stereochimica 3D è ricavata **dalle coordinate**: se sono planari/degenerate, l'app non assegna alcuna configurazione. La visualizzazione e la geometria **non sono** una minimizzazione energetica.

Leggi [Guida alla stereochimica](docs/STEREOCHIMICA_CIP.md) per i casi R/S e [Isomeria E/Z](docs/ISOMERIA_EZ.md) per l’analisi degli alcheni.

## Novità della v0.5: isomeria geometrica E/Z e guida completa

- Individuazione e assegnazione **E/Z** ai doppi legami **C=C in molecole acicliche** con due sostituenti diversi ed espliciti su ciascun carbonio (compresi gli idrogeni).
- Pannello che mostra **priorità CIP 1–2 per entrambe le estremità**, badge C=C: E/Z e etichette direttamente sul modello.
- Comando **Genera l’isomero geometrico opposto**: ricostruisce E ↔ Z mantenendo le lunghezze dei legami; non è una simulazione di rotazione possibile intorno al legame C=C.
- Nuovi esempi (E)/(Z)-but-2-ene, (E)/(Z)-1,2-dicloroetene, (E)/(Z)-acido but-2-enoico.
- **Guida completa e istruzioni integrate**, con otto capitoli espandibili sul lato destro, accessibili anche tramite il pulsante Guida e **disponibili offline** nel file `standalone.html`.
- Test su **180 modelli 3D di riferimento RDKit** (descrittori E/Z assegnati offline), più test di regressione R/S, CJSON, misure e geometria.

**Limiti E/Z:** niente anelli, cumuleni, isotopi, cariche formali, H impliciti, geometrie fortemente distorte e regole CIP avanzate; l'app restituisce *non assegnabile* dove manca affidabilità. Le distanze derivano esclusivamente dalle coordinate, non da una minimizzazione energetica. Leggi [Isomeria E/Z](docs/ISOMERIA_EZ.md) e [Guida completa](docs/GUIDA_COMPLETA.md).

## Funzioni di modellazione e misura

- Costruzione 3D con 118 elementi e legami singoli/doppi/tripli, aggiunta di H semplificata, rotazione, ingrandimento, spostamento e coordinate x/y/z modificabili.
- Misure di **lunghezze di legame (Å), angoli di legame (°) e angoli diedri (°)**; modifica di lunghezze e diedri in catene acicliche con legame centrale singolo.
- Modelli a sfere e bastoncini, sfere illustrative, bastoncini, filiforme; geometria VSEPR euristica per strutture Lewis semplici.
- Importa/esporta **Chemical JSON `.cjson`** per scambiare atomi, coordinate e legami con Avogadro. Importa/esporta anche `.mol`, `.xyz` e progetto `.json` nativo. L'import SDF legge il primo record.
- Annulla/Ripeti, salvataggio locale senza account, nessuna telemetria inserita nell'app.

**Limiti delle misure:** sono valori matematici calcolati sulle coordinate correnti, non risultati sperimentali. Il comando Riordina spaziatura è una regolarizzazione euristica, **non** MMFF/UFF. I formati MOL/XYZ/CJSON importati non preservano automaticamente tutte le proprietà stereochimiche, isotopiche, orbitali o di spettroscopia; conserva i sorgenti originali.

## Prova immediata senza installazione

1. Estrai lo ZIP e apri **`standalone.html`** con Chrome/Edge/Firefox (funziona offline).
2. Scegli **Esempi molecolari → (R)-butan-2-olo**.
3. Osserva il **centro C2: R**. Clicca il centro o il badge `C2: R` per vedere le priorità CIP.
4. Premi **Crea la molecola speculare**: comparirà **C2: S**. Con Annulla tornerai alla molecola precedente.
5. Scegli **(2R,3R)-acido tartarico** per verificare **due centri R**; lo specchio presenta **due centri S**.
6. Per E/Z scegli **(E)-but-2-ene**, clicca il badge **C2=C3: E**, leggi le priorità CIP e premi **Genera l’isomero geometrico opposto**: comparirà **Z**.
7. Per le misure scegli lo strumento **Misura** e seleziona 2, 3 o 4 atomi nell'ordine richiesto.
8. Apri la sezione **Guida completa e istruzioni** nell'app per tutte le procedure senza Internet.

## Sviluppo

Con Node.js 20+:

```bash
npm run check
npm test
npm run build:standalone
```

Per usare la versione modulare `index.html` (richiede un server per i moduli JS):

```bash
python3 -m http.server 8000
```

Quindi apri `http://localhost:8000/`.

I test includono **300 strutture di riferimento R/S e 180 E/Z ottenute indipendentemente con RDKit** in `test/fixtures/cip_rdkit_300.json`, oltre ai test di regressione delle versioni precedenti. Il manuale è `docs/GUIDA_COMPLETA.md`; la documentazione del collaudo manuale è in `docs/TEST_MANUALE.md`.

## Copyright e relazione con Avogadro

Il codice originale è MIT (`LICENSE`). Non sono state incorporate librerie RDKit, Avogadro, Qt né loro eseguibili. Il formato CJSON è implementato indipendentemente su specifiche pubbliche. Il collegamento a [Avogadro](https://avogadro.cc/) indica interoperabilità di file, non affiliazione né approvazione del progetto originale. L'integrazione del suo motore WebAssembly resta una ricerca separata e richiede ulteriori verifiche tecniche e legali (`docs/PIANO_WEBASSEMBLY.md` e `docs/LICENZE_E_LIMITI.md`).

## Distribuzione pubblica e GitHub Pages

Questa è una **copia indipendente della versione 0.5**, predisposta per un repository pubblico senza importare la cronologia Git del repository privato di sviluppo. Il codice e i test di questa copia sono distribuiti con licenza MIT, con i limiti scientifici documentati.

Per pubblicare la web app statica: caricare **i file e le cartelle di questa distribuzione nella radice** di un repository GitHub separato, quindi, solo dopo aver deciso di rendere pubblici repository e sito, aprire **Settings → Pages → Build and deployment → Deploy from a branch → main → /(root) → Save**. Il sito utilizza direttamente `index.html`, `src/` e `.nojekyll`; **non richiede installazione di Node.js, credenziali o servizi di calcolo**. Per l'uso completamente offline aprire `standalone.html`.

I file di questo pacchetto non contengono sorgenti o binari Avogadro/Qt/RDKit; il termine Avogadro indica la sola compatibilità di formati. Il motore stereochimico è didattico e limitato ai casi descritti. La pubblicazione su GitHub Pages non garantisce accessibilità piena né accuratezza per uso professionale.


# MolecolaLab Web 0.5 — guida completa in italiano

Questa guida corrisponde alla sezione **«Guida completa e istruzioni»** integrata nell’app: tutte le indicazioni essenziali rimangono disponibili anche aprendo `standalone.html` senza Internet. Il pulsante **Guida** nella barra superiore apre le stesse istruzioni in una finestra dedicata.

## 1. Avvio e visualizzazione

1. Apri `standalone.html` nel browser (funziona offline). Se sviluppi a partire dai sorgenti, avvia `python3 -m http.server 8000` nella cartella ed entra in `http://localhost:8000`: non aprire `index.html` direttamente come file locale perché i moduli JavaScript potrebbero essere bloccati.
2. Scegli **Esempi molecolari** oppure **Nuovo** per creare una struttura vuota.
3. Usa lo strumento **Ruota**: trascina con il mouse, scorri la rotellina per lo zoom; su touch trascina con un dito e usa due dita per ingrandire. Con doppio clic o **Reimposta la vista** torni all'inquadratura iniziale.
4. Nel menu **Visualizzazione** scegli sfere/bastoncini, filiforme, bastoncini o sfere illustrative.

## 2. Costruire un modellino molecolare

1. Scegli un elemento tra i 118 disponibili, poi lo strumento **Atomo** e clicca lo sfondo.
2. Per creare un nuovo atomo collegato, clicca prima un atomo già presente e poi un punto libero del visore.
3. Per unire due atomi già presenti, scegli **Legame**, imposta singolo/doppio/triplo e clicca due atomi.
4. Con **Seleziona** clicca un atomo per modificarne elemento, carica formale e coordinate x/y/z; clicca un legame per cambiarne ordine o lunghezza.
5. **Sposta** trascina un atomo, **Elimina** rimuove oggetti; **Completa H** e **Riordina spaziatura** sono strumenti euristici, non un calcolo energetico.

## 3. Lunghezze, angoli di legame, diedri

1. Nel pannello **Misure geometriche** scegli la grandezza: **Distanza** (2 atomi, Å), **Angolo** (3 atomi, °) oppure **Diedro** (4 atomi, °).
2. Seleziona **Misura** e clicca gli atomi in sequenza. Nell'angolo A–B–C il vertice è B.
3. Per leggere e cambiare un diedro A–B–C–D, usa **Controllo dell'angolo diedro**. Il legame B–C deve essere singolo, fuori da un anello.
4. I risultati sono geometrici e basati sulle coordinate; non sono previsioni di energia, equilibrio o distanze sperimentali.

## 4. Configurazione R/S

1. Carica **(R)-butan-2-olo**, **(S)-alanina** o uno degli altri esempi R/S.
2. Clicca il centro tetraedrico nel modellino oppure il badge **C…: R/S** nel pannello.
3. Leggi le **priorità CIP 1–4**, poi il descrittore R o S. La rotazione della camera non cambia la configurazione.
4. **Crea la molecola speculare** cambia la disposizione spaziale di tutti i centri stereogenici (R ↔ S per i centri supportati); per tornare indietro usa **Annulla**.

**Ambito R/S:** grafi aciclici, centri C/Si/N⁺ con quattro legami singoli e quattro ligandi espliciti. Non vengono assegnati isotopi, anelli, H impliciti, pseudoasimmetria e regole CIP avanzate.

## 5. Configurazione E/Z (nuovo in 0.5)

1. Scegli dal menu **(E)-but-2-ene**, **(Z)-but-2-ene**, **(E)/(Z)-1,2-dicloroetene** o **(E)/(Z)-acido but-2-enoico**.
2. Clicca con **Ruota** o **Seleziona** al centro del doppio legame C=C; puoi anche cliccare il badge **C2=C3: E/Z** nel pannello **Isomeria geometrica · E/Z**.
3. Leggi le **priorità CIP 1 e 2 su ogni carbonio**. Si confrontano i due sostituenti di priorità 1:
   - **Z** (*zusammen*): stanno dalla **stessa parte** del doppio legame.
   - **E** (*entgegen*): stanno da **parti opposte**.
4. Premi **Genera l'isomero geometrico opposto (E ↔ Z)**. Il software riorienta rigidamente un intero frammento per costruire una struttura dell'altro isomero. **Non simula una rotazione possibile intorno a un doppio legame**, che ha barriera energetica elevata.
5. Usa **Annulla** per ripristinare il modello di partenza. Puoi misurare le distanze prima e dopo: le lunghezze dei legami restano invariate perché la trasformazione è rigida.

**Ambito E/Z:** doppi legami C=C in molecole **acicliche**, su ogni carbonio due ligandi *diversi* con legami singoli ed **esplicitamente disegnati** (inclusi gli H). Coordinate 3D attendibili, non deformate. Sono esclusi anelli, isotopi, gruppi con carica formale, cumuleni, H impliciti e i casi che richiedono regole CIP oltre al confronto dei numeri atomici per livelli. Se i due sostituenti sono equivalenti (come nell'etene), **E/Z non è definito**. Il software non attribuisce E/Z quando non può farlo in modo affidabile.

**Differenza rispetto a R/S:** riflettere l'intera molecola in uno specchio inverte i centri R/S, ma **non** converte E in Z e viceversa. E/Z è una relazione relativa fra le due estremità del doppio legame.

## 6. Previsioni VSEPR e avvertenze sulla geometria

La scheda **Geometria e VSEPR** usa un conteggio elettronico semplificato. Il valore previsto non sostituisce le misure geometriche; dopo ogni modifica controlla esplicitamente le valenze, le cariche e i domini elettronici. Non sono presenti calcoli quantomeccanici, un campo di forza MMFF/UFF in esecuzione, orbitali o predizioni energetiche.

## 7. Salvare e trasferire i file

- **Progetto JSON**: consigliato per salvare integralmente il documento gestito dall'app.
- **Chemical JSON `.cjson`**: scambio di atomi, coordinate e legami di base con Avogadro; alcune proprietà non supportate non vengono mantenute.
- **MOL V2000**: connettività e coordinate di base; non garantisce la conservazione completa delle informazioni stereochimiche.
- **XYZ**: coordinate ed elementi, **senza** la connettività dei legami (che in importazione viene stimata per distanza), quindi non è un archivio stereochimico sicuro.
- **SDF**: l'importazione legge solo il primo record.

Conserva sempre la copia originale dei file importati, soprattutto quando contengono stereochimica o dati scientifici avanzati. Non ci sono invii remoti di molecole incorporati nell'app.

## 8. Scorciatoie, errori e ripristino

- Ctrl/Cmd + Z: annulla; Ctrl/Cmd + Y (o Ctrl/Cmd + Shift + Z): ripeti.
- Esc: deseleziona. Canc/Backspace: elimina l'oggetto selezionato quando non stai scrivendo in un campo.
- Se R/S o E/Z non sono assegnati, leggi il motivo indicato nel pannello; non interpretare `?` come una configurazione valida.
- Se annulli un'inversione E/Z o una modifica della geometria, le coordinate precedenti sono ripristinate dal sistema Annulla/Ripeti.
- Non usare i risultati come certificazione di struttura, identificazione stereochimica completa o calcolo quantitativo senza validazione indipendente.

## Percorsi suggeriti per il liceo scientifico

**Esercizio A (10 minuti):** apri H₂O, misura le due distanze O–H e l'angolo H–O–H; confronta con la previsione VSEPR.

**Esercizio B (15 minuti):** apri (R)-butan-2-olo, annota CIP 1–4 e R/S; crea l'enantiomero e controlla il cambiamento.

**Esercizio C (15 minuti):** confronta (E)- e (Z)-but-2-ene, seleziona C=C, annota CIP 1/2 sulle due estremità; genera l'altro isomero, verifica le lunghezze e spiega perché non si tratta di una rotazione libera del doppio legame.

**Esercizio D (5 minuti):** apri etene C₂H₄, verifica perché **E/Z non è assegnabile** e spiega la necessità di due sostituenti diversi per ciascun carbonio sp².

# Checklist di collaudo — MolecolaLab Web v0.5

## Prove stereochimiche

1. Apri `standalone.html` offline: seleziona (R)-bromoclorofluorometano. Deve comparire `C2: R`. Seleziona C2 e controlla le priorità `Br > Cl > F > H`.
2. Carica la versione (S): il badge deve cambiare in `C2: S`.
3. Carica (R)-butan-2-olo e (S)-butan-2-olo. La prima deve dare R, la seconda S; il C2 ha priorità `OH > etile > metile > H`.
4. Carica (S)-acido lattico: `C2: S` e `OH > COOH > CH3 > H`. Carica (S)-alanina: `C2: S` e `NH2 > COOH > CH3 > H`.
5. Carica acido tartarico (2R,3R): devono comparire **due centri R** (ID 4,6 nell'editor); la funzione speculare deve convertirli entrambi in S.
6. Con (R)-bromoclorofluorometano, premi **Crea la molecola speculare**: `R → S`. Premi Annulla: `S → R`. Ripeti: `R → S`.
7. Ruota la telecamera e nascondi/mostra le etichette. Il descrittore deve rimanere uguale.
8. Esporta (R)-butan-2-olo in JSON nativo e `.cjson`, poi reimporta. Le coordinate e la configurazione R devono rimanere coerenti.
9. Carica metano, benzene e ammonio: non devono apparire centri R/S inesistenti.
10. Importa (se disponibile) una molecola con centro tetraedrico ma H implicito o anello: non deve apparire una configurazione inventata.

## Funzioni precedenti, da preservare

11. H₂O: misura O–H ≈ 0,957 Å e H–O–H ≈ 104,51°. Controlla la previsione VSEPR AX₂E₂.
12. Butano: imposta il diedro degli atomi 1–2–3–4 a 180°, verifica la lettura, annulla/ripeti e controlla che le lunghezze restino costanti.
13. Carbonato: round-trip `.cjson` con ordini 2/1/1 e carica totale −2. Importazioni di CJSON non supportate devono avvisare o bloccare la perdita silenziosa di proprietà.
14. Usa lo schermo stretto e lo zoom touch/pinch in browser emulato; **l'iPad reale non è collaudato**.

## Test automatici

Da terminale nella directory:

```bash
npm run check
npm test
npm run build:standalone
```

La suite comprende i 55 test delle versioni precedenti, 42 test aggiuntivi sulla stereochimica e **300 prove indipendenti su strutture RDKit**, per un totale di 397 test. I dati RDKit sono fixture statiche: non è necessario installare Python o RDKit per eseguire `npm test`. Inoltre, 10 gruppi di prove sono eseguiti con Chromium automatizzato.

La verifica in Chromium non equivale a un test su iPad fisico; non include il motore Avogadro WebAssembly, che non è distribuito nella v0.5.


## Collaudo v0.5 — E/Z e manuale integrato

1. In esempi caricare (E)- e (Z)-but-2-ene: il pannello deve indicare rispettivamente E e Z; selezionando il badge deve mostrare priorità CIP 1 e 2 per C2 e C3.
2. Generare l'isomero geometrico opposto: E diventa Z (e viceversa), le lunghezze di tutti i legami restano identiche; **Annulla** riporta l'isomero iniziale e **Ripeti** lo recupera.
3. Aprire (E)- e (Z)-1,2-dicloroetene e acido but-2-enoico, con visualizzazione e spiegazione coerenti.
4. Con etene C₂H₄ controllare che non venga inventata una configurazione E/Z (i sostituenti su ciascun carbonio sono equivalenti).
5. Il comando **Molecola speculare** deve invertire R/S quando applicabile ma non deve modificare E/Z del but-2-ene.
6. Distorgere le coordinate o cancellare gli H di un esempio E/Z: il pannello deve riportare un motivo di non assegnazione invece di un E/Z fittizio.
7. Premere **Guida**: verificare che la finestra abbia otto capitoli, si possa scorrere, aprire le sezioni e chiudere. La stessa guida deve apparire nel pannello laterale.
8. Verificare il file `standalone.html` offline, l'uso su viewport desktop e mobile, il pinch, e la conservazione E/Z nei salvataggi di progetto JSON e CJSON.

Casi non coperti: atomi isotopici, anelli, cumuleni, H impliciti e stereochimica condizionata da altre unità stereogeniche.

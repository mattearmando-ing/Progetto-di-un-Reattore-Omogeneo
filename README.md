# Progetto di un Reattore Omogeneo per la Produzione di Radionuclidi (Mo-99)

![Status](https://img.shields.io/badge/Status-Completed-brightgreen)
![Domain](https://img.shields.io/badge/Domain-Nuclear%20Engineering-blue)
![Simulation](https://img.shields.io/badge/Tool-MATLAB-orange)

## 📖 Overview
L'obiettivo di questo progetto è lo sviluppo di un modello per un reattore nucleare omogeneo a soluzione acquosa. Il reattore è destinato alla produzione di radionuclidi, in particolare il Molibdeno-99 (Mo-99). Il Mo-99 è di primaria importanza per le applicazioni in ambito medico (medicina nucleare per la produzione di Tecnezio-99m) e industriale. 

Il lavoro analizza le grandezze fisiche e ingegneristiche necessarie per la progettazione e l'esercizio del reattore, includendo il calcolo del rapporto di moderazione ottimale, della sovracriticità iniziale e l'evoluzione temporale della produzione isotopica.

## 🔬 Key Features
- **Combustibile:** Solfato di uranile ($UO_2SO_4$) in soluzione acquosa.
- **Isotopo Fissile:** Uranio-235 arricchito al 19,75% (LEU - Low Enriched Uranium).
- **Moderatore e Refrigerante:** Acqua leggera ($H_2O$).
- **Geometria:** Cilindro equilatero ottimizzato per minimizzare le perdite neutroniche, con volume utile di 30 litri (Raggio = 16,84 cm, Altezza = 33,68 cm).
- **Condizioni Operative:** Potenza termica di progetto pari a 200 kW con una temperatura di funzionamento di 70°C.

## 📊 Key Results
| Metrica | Valore Calcolato |
| :--- | :--- |
| **Rapporto di moderazione ($\theta$)** | `50` (condizione sottomoderata per sicurezza) |
| **Moltiplicazione infinita ($k_{\infty}$)** | `1.05416` |
| **Controllo Reattività ($k_{eff} = 1$)** | `47,78 g/L` di acido borico |
| **Flusso Neutronico (Max)** | `2,2039e13 neutroni/cm²·s` |
| **Vita utile combustibile** | `~81 giorni` (prima che $k_{\infty}$ scenda sotto 1) |

> **⚠️ Note sulla Produzione:** Lo studio mette a confronto la produzione di Mo-99 tramite fissione nucleare dell'U-235 (metodo principale) e tramite attivazione neutronica del Mo-98. La fissione garantisce un tasso di produzione notevolmente superiore, raggiungendo la saturazione in tempi più brevi rispetto all'attivazione.

## 🛠️ Technologies & Tools
- **Analisi e Simulazione:** MATLAB (per il calcolo della formula a quattro fattori, l'evoluzione del combustibile e le curve di decadimento).
- **Librerie Nucleari:** Dati estratti da JEFF-3.3 per le sezioni d'urto.

## 📄 Full Report
Per una descrizione dettagliata del dimensionamento, dell'analisi parametrica e dei modelli di deplezione, consulta il documento di riferimento: **Progetto_reattore_omogeneo.pdf**.

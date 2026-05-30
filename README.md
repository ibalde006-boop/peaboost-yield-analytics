# 🔬 Pipeline de Simulation Variétale & Phénotypage Drone (Pois Protéagineux)

Ce dépôt héberge l'étude de cas technique développée en autonomie dans le cadre du **Master 2 ASAD (USSEIN)**[cite: 1]. Face à la restructuration du stage PeaBoost 2026, ce projet démontre la mise en place d'un pipeline complet d'analyse de données agronomiques et de traitement d'imagerie drone à partir de données simulées biologiquement réalistes.

## 📊 Structure de l'Étude de Cas
Le pipeline R est structuré en 6 grandes étapes scientifiques :
1. **Simulation Expérimentale :** Génération de données synthétiques (8 génotypes x 4 environnements x 2 campagnes) intégrant composantes du rendement (PMG, LAI_max, protéines) et indices spectraux.
2. **Modèles Mixtes (`lme4`) :** Modélisation des effets G×E et classification des moyennes marginales par test de Tukey.
3. **Analyse G×E de Stabilité (`metan`) :** Modèle AMMI, indice WAASB et graphiques *GGE Biplot*.
4. **Analyse Multivariée (`FactoMineR`) :** ACP et Classification Ascendante Hiérarchique (CAH) pour le groupement des profils agronomiques.
5. **Télédétection Drone (`terra` & `RStoolbox`) :** Simulation d'une matrice multibande (Micasense RedEdge), calcul d'indices spectraux (**NDVI, GNDVI, GRVI, NDRE, SAVI**) et régressions linéaires prédictives.
6. **Aide à la Décision :** Tableau de synthèse final pour la sélection variétale.

## 🛠️ Stack R utilisée
`lme4` | `metan` | `FactoMineR` | `terra` | `RStoolbox` | `tidyverse` | `ggplot2` | `corrplot`

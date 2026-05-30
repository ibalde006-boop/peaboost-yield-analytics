# ==============================================================================
# ÉTUDE DE CAS COMPLÈTE - SIMULATION ET ANALYSE VARIÉTALE PEABOOST 2026
# Analyse de l'élaboration du rendement du pois protéagineux
# Auteur: Ibrahima BALDE (Master 2 ASAD - USSEIN)
# ==============================================================================

# ÉTAPE 0 - Chargement des packages principaux
library(lme4)       # Modèles mixtes
library(lmerTest)   # Tests de significativité des effets fixes
library(emmeans)    # Comparaisons de moyennes et CLD
library(metan)      # Analyse de l'interaction GxE, AMMI et GGE biplot
library(FactoMineR) # ACP et Classification (HCPC)
library(factoextra) # Visualisation élégante de l'ACP
library(terra)      # Manipulation d'images matricielles (Drone/Sat)
library(RStoolbox)  # Calculs avancés d'indices spectraux
library(tidyverse)  # Data manipulation (dplyr, tidyr...)
library(ggplot2)    # Dataviz avancée
library(corrplot)   # Matrices de corrélations graphiques

# ==============================================================================
# PARTIE 1 - SIMULATION DES DONNÉES EXPÉRIMENTALES
# ==============================================================================
set.seed(42)

# Définition des facteurs (8 génotypes x 4 stations x 3 blocs x 2 ans)
genotypes <- c("G1_Kayanne", "G2_Isard", "G3_Enduro", "G4_Astronaute", 
               "G5_Baccara", "G6_Champagne", "G7_Lumina", "G8_Diver")
environnements <- c("E1_Toulouse", "E2_Vendôme", "E3_Dijon", "E4_Evreux")

donnees_pois <- expand.grid(
  Genotype = genotypes,
  Environnement = environnements,
  Repetition = 1:3,
  Campagne = c(2025, 2026)
)

# Définition des potentiels biologiques (effets théoriques)
effet_genotype <- c(4.2, 3.8, 4.5, 4.0, 3.5, 4.8, 3.9, 4.3)
names(effet_genotype) <- genotypes
effet_env <- c(0.5, -0.3, 0.2, -0.1)
names(effet_env) <- environnements

# Construction mathématique de la matrice agronomique
donnees_pois <- donnees_pois %>%
  mutate(
    # Rendement grain simulé avec résidu gaussien (t/ha)
    Rendement = effet_genotype[Genotype] + effet_env[Environnement] + rnorm(n(), 0, 0.3),
    # Composantes du rendement liées de manière logique au rendement global
    Nb_gousses_m2 = round(80 + Rendement * 12 + rnorm(n(), 0, 8)),
    Nb_graines_gousse = round(3.2 + rnorm(n(), 0, 0.4), 1),
    PMG = round(200 + Rendement * 8 + rnorm(n(), 0, 15)),
    # Variables phénologiques (Jours Après Levée)
    Floraison_JAL = round(110 + rnorm(n(), 0, 5)),
    Maturite_JAL = round(170 + rnorm(n(), 0, 7)),
    # Indicateurs physiologiques et spectraux
    LAI_max = round(3.5 + Rendement * 0.2 + rnorm(n(), 0, 0.3), 2),
    Proteines = round(22 + rnorm(n(), 0, 1.5), 1),
    NDVI_floraison = round(0.65 + Rendement * 0.03 + rnorm(n(), 0, 0.04), 3),
    GNDVI_floraison = round(0.55 + Rendement * 0.025 + rnorm(n(), 0, 0.035), 3)
  )

# ==============================================================================
# PARTIE 2 - MODÈLES MIXTES AVEC lme4
# ==============================================================================
cat("\n", rep("=", 60), "\n")
cat("PARTIE 2 - MODÈLES MIXTES AVEC lme4\n")
cat(rep("=", 60), "\n\n")

modele_mixte <- lmer(
  Rendement ~ Genotype + Environnement + Campagne + Genotype:Environnement + (1 | Environnement:Repetition),
  data = donnees_pois
)

cat(" - Résumé du modèle mixte -\n")
print(summary(modele_mixte))

cat("\n - ANOVA Type III (lmerTest) -\n")
print(anova(modele_mixte))

cat("\n - Moyennes marginales estimées par génotype -\n")
emm_genotype <- emmeans(modele_mixte, ~ Genotype)
print(emm_genotype)

cat("\n - Test de comparaison multiple (Tukey) -\n")
print(pairs(emm_genotype, adjust = "tukey"))

# ==============================================================================
# PARTIE 3 - ANALYSE GxE AVEC metan
# ==============================================================================
cat("\n", rep("=", 60), "\n")
cat("PARTIE 3 - ANALYSE GxE AVEC metan\n")
cat(rep("=", 60), "\n\n")

donnees_gxe <- donnees_pois %>%
  group_by(Genotype, Environnement) %>%
  summarise(Rendement = mean(Rendement), .groups = "drop")

cat(" - Modèle AMMI -\n")
modele_ammi <- performs_ammi(donnees_gxe, env = Environnement, gen = Genotype, rep = NULL, resp = Rendement)
print(modele_ammi)

cat("\n - Indice de Stabilité des génotypes (WAASB) -\n")
stabilite <- waasb(donnees_gxe, env = Environnement, gen = Genotype, rep = NULL, resp = Rendement)
print(stabilite)

cat("\n - Modélisation GGE Biplot -\n")
gge_model <- gge(donnees_gxe, env = Environnement, gen = Genotype, resp = Rendement)
print(gge_model)

# ==============================================================================
# PARTIE 4 - ACP ET TYPOLOGIES AVEC FactoMineR
# ==============================================================================
cat("\n", rep("=", 60), "\n")
cat("PARTIE 4 - ACP ET TYPOLOGIES AVEC FactoMineR\n")
cat(rep("=", 60), "\n\n")

donnees_acp <- donnees_pois %>%
  group_by(Genotype) %>%
  summarise(
    Rendement = mean(Rendement), Nb_gousses_m2 = mean(Nb_gousses_m2),
    Nb_graines_gousse = mean(Nb_graines_gousse), PMG = mean(PMG),
    Floraison_JAL = mean(Floraison_JAL), Maturite_JAL = mean(Maturite_JAL),
    LAI_max = mean(LAI_max), Proteines = mean(Proteines),
    NDVI_floraison = mean(NDVI_floraison), GNDVI_floraison = mean(GNDVI_floraison),
    .groups = "drop"
  )

mat_acp <- donnees_acp %>% column_to_rownames("Genotype") %>% as.data.frame()

res_acp <- PCA(mat_acp, scale.unit = TRUE, graph = FALSE)
cat(" - Variance expliquée par axe -\n")
print(res_acp$eig)

cat("\n - Contributions des variables aux axes -\n")
print(res_acp$var$contrib)

cat("\n - Classification Ascendante Hiérarchique (CAH) -\n")
res_hcpc <- HCPC(res_acp, nb.clust = 3, graph = FALSE)
cat("\n - Description des clusters par les variables -\n")
print(res_hcpc$desc.var)

# ==============================================================================
# PARTIE 5 - INDICES SPECTRAUX DRONE AVEC terra
# ==============================================================================
cat("\n", rep("=", 60), "\n")
cat("PARTIE 5 - INDICES SPECTRAUX DRONE AVEC terra\n")
cat(rep("=", 60), "\n\n")

# Simulation d'une matrice d'imagerie multibande drone (50x50 px, 5 bandes)
image_drone <- rast(nrows = 50, ncols = 50, nlyr = 5)
values(image_drone[[1]]) <- runif(2500, 0.05, 0.15) # Bleu
values(image_drone[[2]]) <- runif(2500, 0.10, 0.25) # Vert
values(image_drone[[3]]) <- runif(2500, 0.08, 0.20) # Rouge
values(image_drone[[4]]) <- runif(2500, 0.20, 0.45) # RedEdge
values(image_drone[[5]]) <- runif(2500, 0.40, 0.75) # PIR
names(image_drone) <- c("Bleu", "Vert", "Rouge", "RedEdge", "PIR")

# Calcul programmatique des équations d'indices de réflectance
NDVI <- (image_drone[["PIR"]] - image_drone[["Rouge"]]) / (image_drone[["PIR"]] + image_drone[["Rouge"]])
GNDVI <- (image_drone[["PIR"]] - image_drone[["Vert"]]) / (image_drone[["PIR"]] + image_drone[["Vert"]])
GRVI <- (image_drone[["Vert"]] - image_drone[["Rouge"]]) / (image_drone[["Vert"]] + image_drone[["Rouge"]])
NDRE <- (image_drone[["PIR"]] - image_drone[["RedEdge"]]) / (image_drone[["PIR"]] + image_drone[["RedEdge"]])
L <- 0.5
SAVI <- ((image_drone[["PIR"]] - image_drone[["Rouge"]]) / (image_drone[["PIR"]] + image_drone[["Rouge"]] + L)) * (1 + L)

indices_stack <- c(NDVI, GNDVI, GRVI, NDRE, SAVI)
names(indices_stack) <- c("NDVI", "GNDVI", "GRVI", "NDRE", "SAVI")

cat(" - Statistiques globales des rasters calculés -\n")
print(global(indices_stack, fun = c("mean", "sd")))

# Modélisation / Régression linéaire
correlation_data <- donnees_pois %>%
  group_by(Genotype) %>%
  summarise(Rendement = mean(Rendement), NDVI_drone = mean(NDVI_floraison), .groups = "drop")

modele_ndvi <- lm(Rendement ~ NDVI_drone, data = correlation_data)
cat("\n - Analyse de la régression linéaire Rendement ~ NDVI Drone -\n")
print(summary(modele_ndvi))

# ==============================================================================
# PARTIE 6 - SYNTHÈSE ET AIDE À LA DÉCISION
# ==============================================================================
cat("\n", rep("=", 60), "\n")
cat("PARTIE 6 - SYNTHÈSE et AIDE À LA DÉCISION SÉLECTION\n")
cat(rep("=", 60), "\n\n")

synthese_finale <- donnees_pois %>%
  group_by(Genotype) %>%
  summarise(
    Rendement_moyen = round(mean(Rendement), 2),
    PMG_moyen = round(mean(PMG), 0),
    Nb_gousses_moyen = round(mean(Nb_gousses_m2), 0),
    NDVI_moyen = round(mean(NDVI_floraison), 3),
    Proteines_moy = round(mean(Proteines), 1),
    LAI_max_moy = round(mean(LAI_max), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(Rendement_moyen))

print(as.data.frame(synthese_finale))

cat("\n[SUCCÈS] Traitement terminé. Variété élite identifiée :", synthese_finale$Genotype[1], "\n")

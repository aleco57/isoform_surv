### PCA and UMAP ###
load("Results/ciso_pheno_surv.RData")
iso_exp <- BRCA_surv[,-c(1:10)]

## PCA
library(tidymodels)

pca_rec <- recipe(~., data = iso_exp) %>% 
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors())
pca_prep <- prep(pca_rec)
pca_prep

tidied_pca <- tidy(pca_prep, 2)

#Overview of first 5 PCs
tidied_pca %>%
  filter(component %in% paste0("PC", 1:5)) %>%
  mutate(component = forcats::fct_inorder(component)) %>%
  ggplot(aes(value, terms, fill = terms)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~component, nrow = 1) +
  labs(y = NULL)

#We can zoom in to look at those 8 isoforms which contribute most to the first four PCA
library(tidytext)
tidied_pca %>%
  filter(component %in% paste0("PC", 1:4)) %>%
  group_by(component) %>%
  top_n(8, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(
    x = "Absolute value of contribution",
    y = NULL, fill = "Positive?"
  )


juice(pca_prep)


#OR WE COULD USE ANOTHER METHOD FOR PCA:
pc <- prcomp(BRCA_surv[,-c(1:10)], scale = TRUE)
summary(pc)
pc$rotation

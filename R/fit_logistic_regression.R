logistic_regression <- function(x_train, y_train, x_test, y_test, perform_lasso = FALSE) {
  library(glmnet)
  library(caret)
  library(pROC)
  x_train <- scale(t(log2(x_train+0.01)))
  x_test <- scale(t(log2(x_test+0.01)))
  x_train <- x_train[ , colSums(is.na(x_train))==0]
  x_test <- x_test[ , colSums(is.na(x_test))==0]
  common_features <- intersect(colnames(x_train), colnames(x_test))
  x_train <- x_train[, common_features]
  x_test <- x_test[, common_features]
  x_train_matrix <- as.matrix(x_train)
  x_test_matrix <- as.matrix(x_test)
  y_train <- ifelse(y_train == "group_1", 1, 0)
  y_test <- ifelse(y_test == "group_1", 1, 0)
  if (perform_lasso) {
    cv.fit = cv.glmnet(x_train_matrix, y_train, family = "binomial", type.measure = "class", alpha = 1)
    best_lambda = cv.fit$lambda.min
    fit = glmnet(x_train_matrix, y_train, family = "binomial", alpha = 1, lambda = best_lambda)
    preds = predict(fit, x_test_matrix, type = "response")[,1]
  } else {
    fit = glm(y_train ~ ., family = binomial(), data = as.data.frame(x_train))
    preds = predict(fit, as.data.frame(x_test), type = "response")
  }
  pred_classes = ifelse(preds > 0.5, 1, 0)
  predictions_df <- data.frame(preds, pred_classes, y_test)
  colnames(predictions_df) <- c("pred_probs", "pred_labels", "true_labels")
  y_levels <- levels(factor(y_test))
  y_pred_factor <- factor(pred_classes, levels = y_levels)
  y_test_factor <- factor(y_test, levels = y_levels)
  performance_metrics <- confusionMatrix(y_pred_factor, y_test_factor)
  roc_obj <- roc(y_test, pred_classes)
  roc_auc <- auc(roc_obj)
  metrics <- list(sensitivity = performance_metrics$byClass[1],
                  specificity = performance_metrics$byClass[2],
                  balanced_accuracy = performance_metrics$byClass[11],
                  auc = roc_auc)
  coef_df <- as.data.frame(as.matrix(coef(fit)))
  return(list(logistic_metrics = metrics, logistic_coef = coef_df, logistic_predictions = predictions_df))
}

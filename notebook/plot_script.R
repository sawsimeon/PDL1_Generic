source("script.R")

get_recall_precision = function(precision_recall) {
  data_df = precision_recall$curve
  recall = data_df[, 1]
  precision = data_df[, 2]
  result_df = data.frame(precision = precision,
                         recall = recall)
  return(result_df)
}
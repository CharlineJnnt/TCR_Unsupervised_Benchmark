# similarity functions

jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  union = length(A) + length(B) - intersection
  return (intersection/union)
}

perform_jaccard_cluster <- function(list_alignement_1, list_alignement_2){
  edge_cluster_jaccard <- data.frame()
  for(cluster_1 in names(list_alignement_1)){
    for(cluster_2 in names(list_alignement_2)){
      score <- jaccard_similarity(list_alignement_1[[cluster_1]]$CDR3, list_alignement_2[[cluster_2]]$CDR3)
      size_1 <- nrow(list_alignement_1[[cluster_1]])
      size_2 <- nrow(list_alignement_2[[cluster_2]])
      df <- data.frame(cluster_1 = cluster_1, cluster_2 = cluster_2, jaccard_index = score, size_1 = size_1, size_2 = size_2)
      edge_cluster_jaccard <- rbind(edge_cluster_jaccard, df)
    }
  }
  return(edge_cluster_jaccard)
}

# overlap_coefficient <- function(A, B) {
#   intersection = length(intersect(A, B))
#   min = min(length(A), length(B))
#   return (intersection/min)
# }
# 
# perform_overlap_coeff_cluster <- function(list_alignement_1, list_alignement_2){
#   edge_cluster_overlap <- data.frame()
#   for(cluster_1 in names(list_alignement_1)){
#     for(cluster_2 in names(list_alignement_2)){
#       score <- overlap_coefficient(list_alignement_1[[cluster_1]]$CDR3, list_alignement_2[[cluster_2]]$CDR3)
#       size_1 <- nrow(list_alignement_1[[cluster_1]])
#       size_2 <- nrow(list_alignement_2[[cluster_2]])
#       df <- data.frame(cluster_1 = cluster_1, cluster_2 = cluster_2, overlap_index = score, size_1 = size_1, size_2 = size_2)
#       edge_cluster_overlap <- rbind(edge_cluster_overlap, df)
#     }
#   }
#   return(edge_cluster_overlap)
# }
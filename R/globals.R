# Define global variables to satisfy R CMD check for ggplot2 aesthetics and internal bindings
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("is_admissible", "kD", "kY", "admissible_num", "dof1"))
}

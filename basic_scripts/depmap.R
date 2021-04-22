library("depmap")
library("ExperimentHub")

## create ExperimentHub query object
eh <- ExperimentHub()
query(eh, "depmap")

## retrieve metadata about cancer cell lines
metadata <- depmap::depmap_metadata()

print(metadata)

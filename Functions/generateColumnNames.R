### generateColumnNames(plateName)
### 
### Parameters
###   - plateName: Pass a string plateName, example: "JE1"
### 
### Return
###   - 
### 
### Authors
###   - Joep Eding
### 
### TODO
###   - 
### Returns a list (dataframe?) of column names generated as `platename`_`row letter`_`column number` that can be assigned to columnames(dataframe)
generateColumnNames <- function(plateName){
  columnNames <- sapply(c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"), {function(x) sapply(1:24,function(y) {paste(plateName,x,y, sep="_")})})
  return(columnNames)
}
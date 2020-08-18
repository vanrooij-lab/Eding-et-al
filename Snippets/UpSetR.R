if("UpSetR" %in% rownames(installed.packages()) == FALSE) {install.packages("UpSetR")}
library("UpSetR")

topXGenes=100

listInput <- list(
  all.3 = (clusterDifferentialExpression$patientAll$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes))$rowname,
  p1.1 = (clusterDifferentialExpression$patient1$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname,
  p2.4 = (clusterDifferentialExpression$patient2$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname,
#  p4.2 = (clusterDifferentialExpression$patient4$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname,
  p4.6 = (clusterDifferentialExpression$patient4$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname
#  p1.2 = (clusterDifferentialExpression$patient1$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname
)
listInput <- list(
  all.7 = (clusterDifferentialExpression$patientAll$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.5 = (clusterDifferentialExpression$patient1$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname,
  #p2.9 = (clusterDifferentialExpression$patient2$cl.9 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname,
  p2.10 = (clusterDifferentialExpression$patient2$cl.10 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname,
  p4.11 = (clusterDifferentialExpression$patient4$cl.11 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc) %>% select(rowname))$rowname
)
listInput <- list(
  pool.1 = (clusterDifferentialExpression$patientAll$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.2 = (clusterDifferentialExpression$patientAll$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.3 = (clusterDifferentialExpression$patientAll$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.4 = (clusterDifferentialExpression$patientAll$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.5 = (clusterDifferentialExpression$patientAll$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.6 = (clusterDifferentialExpression$patientAll$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.7 = (clusterDifferentialExpression$patientAll$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.8 = (clusterDifferentialExpression$patientAll$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.1 = (clusterDifferentialExpression$patient1$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.2 = (clusterDifferentialExpression$patient1$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.3 = (clusterDifferentialExpression$patient1$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.4 = (clusterDifferentialExpression$patient1$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.5 = (clusterDifferentialExpression$patient1$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.6 = (clusterDifferentialExpression$patient1$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.7 = (clusterDifferentialExpression$patient1$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.8 = (clusterDifferentialExpression$patient1$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.1 = (clusterDifferentialExpression$patient2$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.2 = (clusterDifferentialExpression$patient2$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.3 = (clusterDifferentialExpression$patient2$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.4 = (clusterDifferentialExpression$patient2$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.5 = (clusterDifferentialExpression$patient2$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.6 = (clusterDifferentialExpression$patient2$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.7 = (clusterDifferentialExpression$patient2$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.8 = (clusterDifferentialExpression$patient2$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.9 = (clusterDifferentialExpression$patient2$cl.9 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.10 = (clusterDifferentialExpression$patient2$cl.10 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.1 = (clusterDifferentialExpression$patient4$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.2 = (clusterDifferentialExpression$patient4$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.3 = (clusterDifferentialExpression$patient4$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.4 = (clusterDifferentialExpression$patient4$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.5 = (clusterDifferentialExpression$patient4$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.6 = (clusterDifferentialExpression$patient4$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.7 = (clusterDifferentialExpression$patient4$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.8 = (clusterDifferentialExpression$patient4$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.9 = (clusterDifferentialExpression$patient4$cl.9 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.10 = (clusterDifferentialExpression$patient4$cl.10 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname
)
upset(
  fromList(listInput[seq(length(listInput),1)]), 
#  intersections = list(
#    list('all.1', 'p1.6'),
#    list('all.3', 'p1.1')
#  ),
  order.by = "freq", 
  empty.intersections = NULL, 
  group.by = 'sets', 
  nsets = 36, 
  nintersects = 108,
  mb.ratio = c(0.5,0.5),
  cutoff = 3
)

pdf('UpSet1.pdf')
listInput <- list(
  pool.1 = (clusterDifferentialExpression$patientAll$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.2 = (clusterDifferentialExpression$patientAll$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.3 = (clusterDifferentialExpression$patientAll$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.4 = (clusterDifferentialExpression$patientAll$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.5 = (clusterDifferentialExpression$patientAll$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.6 = (clusterDifferentialExpression$patientAll$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.7 = (clusterDifferentialExpression$patientAll$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.8 = (clusterDifferentialExpression$patientAll$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname
)
UpSetData <- fromList(listInput[seq(length(listInput),1)])
UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  UpSetData,
  empty.intersections = NULL, 
  order.by = "freq", 
  group.by = "degree", 
  nsets = 36, 
  nintersects = 108,
  mb.ratio = c(0.5,0.5),
  queries = list(
    list(query = intersects, params = list('pool.3'), color = "orange", active = T)
  )
)
dev.off()

pdf('UpSet2.pdf')
listInput <- list(
  p1.1 = (clusterDifferentialExpression$patient1$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.4 = (clusterDifferentialExpression$patient2$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.6 = (clusterDifferentialExpression$patient4$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname
)
UpSetData <- fromList(listInput[seq(length(listInput),1)])
UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  UpSetData,
  empty.intersections = NULL, 
  order.by = "freq", 
  group.by = "degree", 
  nsets = 36, 
  nintersects = 108,
  mb.ratio = c(0.5,0.5),
  queries = list(
    list(query = intersects, params = list('p1.1','p2.4','p4.6'), color = "green", active = T)
  )
)
dev.off()

pdf('UpSet3.pdf')
listInput <- list(
  pool.3 = (clusterDifferentialExpression$patientAll$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.1 = (clusterDifferentialExpression$patient1$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.4 = (clusterDifferentialExpression$patient2$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.6 = (clusterDifferentialExpression$patient4$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname
)
UpSetData <- fromList(listInput[seq(length(listInput),1)])
UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  UpSetData,
  empty.intersections = NULL, 
  order.by = "freq", 
  group.by = "degree", 
  nsets = 36, 
  nintersects = 108,
  mb.ratio = c(0.5,0.5),
  queries = list(
    list(query = intersects, params = list('pool.3','p1.1','p2.4','p4.6'), color = "green", active = T),
    list(query = intersects, params = list('pool.3'), color = "orange", active = T)
  )
)
dev.off()

pdf('UpSet4.pdf')
listInput <- list(
  pool.1 = (clusterDifferentialExpression$patientAll$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.2 = (clusterDifferentialExpression$patientAll$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.3 = (clusterDifferentialExpression$patientAll$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.4 = (clusterDifferentialExpression$patientAll$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.5 = (clusterDifferentialExpression$patientAll$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.6 = (clusterDifferentialExpression$patientAll$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.7 = (clusterDifferentialExpression$patientAll$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  pool.8 = (clusterDifferentialExpression$patientAll$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.1 = (clusterDifferentialExpression$patient1$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.2 = (clusterDifferentialExpression$patient1$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.3 = (clusterDifferentialExpression$patient1$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.4 = (clusterDifferentialExpression$patient1$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.5 = (clusterDifferentialExpression$patient1$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.6 = (clusterDifferentialExpression$patient1$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.7 = (clusterDifferentialExpression$patient1$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.8 = (clusterDifferentialExpression$patient1$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.1 = (clusterDifferentialExpression$patient2$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.2 = (clusterDifferentialExpression$patient2$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.3 = (clusterDifferentialExpression$patient2$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.4 = (clusterDifferentialExpression$patient2$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.5 = (clusterDifferentialExpression$patient2$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.6 = (clusterDifferentialExpression$patient2$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.7 = (clusterDifferentialExpression$patient2$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.8 = (clusterDifferentialExpression$patient2$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.9 = (clusterDifferentialExpression$patient2$cl.9 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.10 = (clusterDifferentialExpression$patient2$cl.10 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.1 = (clusterDifferentialExpression$patient4$cl.1 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.2 = (clusterDifferentialExpression$patient4$cl.2 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.3 = (clusterDifferentialExpression$patient4$cl.3 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.4 = (clusterDifferentialExpression$patient4$cl.4 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.5 = (clusterDifferentialExpression$patient4$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.6 = (clusterDifferentialExpression$patient4$cl.6 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.7 = (clusterDifferentialExpression$patient4$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.8 = (clusterDifferentialExpression$patient4$cl.8 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.9 = (clusterDifferentialExpression$patient4$cl.9 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.10 = (clusterDifferentialExpression$patient4$cl.10 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname
)
UpSetData <- fromList(listInput[seq(length(listInput),1)])
UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  UpSetData,
  empty.intersections = NULL, 
  order.by = "freq", 
  group.by = "degree", 
  nsets = 36, 
  nintersects = 252,
  mb.ratio = c(0.5,0.5),
  queries = list(
    list(query = intersects, params = list('pool.8','p1.8','p2.1','p4.9'), color = "red", active = T),
    list(query = intersects, params = list('pool.8','p2.1','p4.9'), color = "red", active = T),
    list(query = intersects, params = list('pool.8','p1.8','p2.1'), color = "red", active = T),
    list(query = intersects, params = list('pool.8','p2.1'), color = "red", active = T),
    list(query = intersects, params = list('pool.8','p1.8'), color = "red", active = T),
    list(query = intersects, params = list('pool.7','p1.5','p2.9','p2.10'), color = "blue", active = T),
    list(query = intersects, params = list('pool.7','p1.5','p2.10'), color = "blue", active = T),
    list(query = intersects, params = list('pool.7','p2.10'), color = "blue", active = T),
    list(query = intersects, params = list('pool.3','p1.1','p2.4','p4.6'), color = "green", active = T),
    list(query = intersects, params = list('pool.3','p4.6'), color = "green", active = T),
    #list(query = intersects, params = list('pool.3','p2.4'), color = "green", active = T),
    list(query = intersects, params = list('pool.3','p1.1'), color = "green", active = T),
    list(query = intersects, params = list('pool.3'), color = "orange", active = T)
  )
)
dev.off()

pdf('KCNQ1OT1_topX.pdf')
listInput <- list(
  #pool.7 = (clusterDifferentialExpression$patientAll$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.5 = (clusterDifferentialExpression$patient1$cl.5 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p2.10 = (clusterDifferentialExpression$patient2$cl.10 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p4.11 = (clusterDifferentialExpression$patient4$cl.11 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname
)
UpSetData <- fromList(listInput[seq(length(listInput),1)])
UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  UpSetData,
  empty.intersections = NULL, 
  order.by = "freq", 
  group.by = "degree", 
  nsets = 36, 
  nintersects = 252,
  mb.ratio = c(0.5,0.5),
  queries = list(
    list(query = intersects, params = list('p1.5','p2.10','p4.11'), color = "green", active = T)
  )
)
dev.off()


pdf('KCNQ1OT1_1.5fc.pdf')
listInput <- list(
  #pool.7 = (clusterDifferentialExpression$patientAll$cl.7 %>% tibble::rownames_to_column() %>% top_n(topXGenes, fc))$rowname,
  p1.5 = (clusterDifferentialExpression$patient1$cl.5 %>% tibble::rownames_to_column() %>% filter(fc > 3))$rowname,
  p2.10 = (clusterDifferentialExpression$patient2$cl.10 %>% tibble::rownames_to_column() %>% filter(fc > 3))$rowname,
  p4.11 = (clusterDifferentialExpression$patient4$cl.11 %>% tibble::rownames_to_column() %>% filter(fc > 3))$rowname
)
UpSetData <- fromList(listInput[seq(length(listInput),1)])
UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  UpSetData,
  empty.intersections = NULL, 
  order.by = "freq", 
  group.by = "degree", 
  nsets = 36, 
  nintersects = 252,
  mb.ratio = c(0.5,0.5),
  queries = list(
    list(query = intersects, params = list('p1.5','p2.10','p4.11'), color = "green", active = T)
  )
)
dev.off()


#' Process Experiment Description Excel Sheet
process.desc <- function() {
  desc <- read.csv('dat/SingleCellExperiments_Description.csv')
  # only take the portion we need
  desc <- desc[1:308,]
  # remove rows that don't have experiment info
  desc <- desc[grepl('^\\d{2,3}[A-Z]{1}', desc$`Exp...`),]
  desc <- desc[,c(1, 3:12)]
  names(desc)[1] <- 'Exp'
  rownames(desc) <- NULL
  desc <- melt(desc, id.vars=c('Exp'), variable.name='Channel', value.name='Sample',
               factorsAsStrings=F)
  desc$Exps <- as.character(desc$Exp)
  # channel ID
  desc$ch <- as.numeric(desc$Channel)
  
  desc$Sample <- as.character(desc$Sample)
  desc$Sample[desc$Sample %in% c('', '0', 'Empty', 'N/A', 'PBS')] <- NA
  
  ## Quantity - how many cells were in each channel
  desc$Quantity <- str_extract(desc$Sample, '\\d+(e\\d)?(\\.\\dk?)?')
  # convert scientific notation
  desc$Quantity[grep('\\de\\d', desc$Quantity)] <- 
    as.numeric(desc$Quantity[grep('\\de\\d', desc$Quantity)])
  # convert "k" notation
  desc$Quantity[grep('\\d\\.\\dk', desc$Quantity)] <- 
    as.numeric(str_extract(desc$Quantity[grep('\\d+\\.\\dk', desc$Quantity)], '\\d+\\.\\d')) * 1000
  
  ## Type - what kind of cell it was
  # J = Jurkat
  # H = Hek293
  # U = U937
  # ES|EB = Mouse
  desc$Type <- str_extract(desc$Sample, '[JHUjhu]{1}|([Ee]{1}[SsBb]{1})')
  # Mixed type - more than one type of cell
  # most likely a carrier channel
  desc$Type[grepl('\\&|and', desc$Sample)] <- 'Mixed'
  
  # Uppercase = intact, single cell
  # Lowercase = diluted, lysed mixture
  desc$Diluted <- F
  desc$Diluted[grepl('[jhu]|es|eb', desc$Type)] <- T
  desc$Diluted[is.na(desc$Type)] <- NA
  
  # Re-normalize cell type for easier searching
  desc$Type <- toupper(desc$Type)
  
  return(desc)
}
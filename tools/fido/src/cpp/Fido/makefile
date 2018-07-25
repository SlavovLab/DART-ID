CPPFLAGS=-Wall -ffast-math -pipe -O4 -g
BINPATH=.
VPATH=../../bin
TARGETS=Fido FidoChooseParameters
#TARGETS=ROC ROC-Interpolate ROC50 ROCDifferenceInterpolate WTest MeanStdDev Plottify StringCounter Fido EM Pivdo Fasta2MissingPeptides MCMC2Pivdo CrossValSplitter Rlist2FPLists FPListIntersection FDR FDRDivergence Rank RelativeOptimize SensitivityAtFDR FPlottify FidoDetectability KDTest FidoKnownDetectability FidoMAP FidoChooseParameters
#TARGETS=FS-LP Grid-LP ROC FP-Princeps ROC-Interpolate Simplex GenerateLP GenerateSparseLP GridSimplex ProtProph2ROC ROC50 ROCDifferenceInterpolate AS-LP WTest exe MeanStdDev Plottify StringCounter Fido EM

all: $(TARGETS)

depend:
	g++ -MM *.cpp *.h

clean:
	rm *.o;

KDTest: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainKDMethod.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

exe: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o Random.o main.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

MCMC2Pivdo: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o  GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o PivdoSplitter.o mainMCMC2Pivdo.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

CrossValSplitter: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o  GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o PivdoSplitter.o mainCrossValSplit.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Plottify: Set.o Vector.o Numerical.o Matrix.o Random.o plottify.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FPlottify: Set.o Vector.o Numerical.o Matrix.o Random.o fplottify.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

MSE: Set.o Vector.o Numerical.o Matrix.o Random.o MSE.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FDRDivergence: Set.o Vector.o Numerical.o Matrix.o Random.o fdrDivergence.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

ROCDifferenceInterpolate: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o Random.o ROCInterpolator.o rocDifferenceInterpolate.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

WTest: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o Random.o w_test.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

MeanStdDev: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o Random.o meanStdDev.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FS-LP: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o mainForwardSampler.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

AS-LP: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o AffineScaling.o Random.o mainAffineScaling.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Simplex: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o Simplex.o Random.o mainSimplex.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

GridSimplex: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o Simplex.o GridSimplex.o Random.o mainGridSimplex.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

GenerateLP: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o Random.o generateLinProg.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

GenerateSparseLP: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o Random.o generateSparseLinProg.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Grid-LP: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o GridLP.o mainGridLP.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Rank: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o rank.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

ROC: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o rlistROC.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

SensitivityAtFDR: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o rlistSensitivityAtFDR.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

RelativeOptimize: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o relativeOptimize.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FDR: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o rlistFDR.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Rlist2FPLists: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o rlistFPLists.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FPListIntersection: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o FPListIntersection.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

StringCounter: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o Random.o stringCounter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FP-Princeps: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o firstFPs.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

ROC50: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o roc50.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

ProtProph2ROC: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o GridLP.o mainROC_ProtProph.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

ROC-Interpolate: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o GridLP.o ROCInterpolator.o interpolateROC.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Fasta2MissingPeptides: Set.o Vector.o Numerical.o Matrix.o Random.o FastaReader.o Proteomics.o testFasta.o ProteinBigraph.o BasicBigraph.o BasicGroupBigraph.o   PivdoBigraph.o ProteinIdentifier.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Pivdo: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinBigraph.o BasicBigraph.o BasicGroupBigraph.o   PivdoBigraph.o ProteinIdentifier.o mainPivdo.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Fido: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o  GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFido.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FidoChooseParameters: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o  GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFido_ChooseParameters.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FidoMAP: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o  GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFidoMAP.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

Fido_Likelihood: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o  GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFido_Likelihood.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FidoKnownDetectability: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o  GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFidoLoadDetectability.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FidoDetectability: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFidoDetectability.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

EM: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o EMBigraph.o  mainEM.o 
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

%.o: %.cpp
	g++ $(CPPFLAGS) -c $<

#automade
AffineScaling.o: AffineScaling.cpp AffineScaling.h LPSolver.h Numerical.h \
  GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp \
  Random.h
Array.o: Array.cpp
BasicBigraph.o: BasicBigraph.cpp BasicBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
BasicGroupBigraph.o: BasicGroupBigraph.cpp BasicGroupBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h
EMBigraph.o: EMBigraph.cpp EMBigraph.h BasicGroupBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h
EmissionModel.o: EmissionModel.cpp EmissionModel.h Combinatorics.h
FastaReader.o: FastaReader.cpp FastaReader.h Array.h Streamable.h Set.h \
  Array.cpp
fdrDivergence.o: fdrDivergence.cpp Matrix.h Vector.h Array.h Streamable.h \
  Set.h Array.cpp Numerical.h
firstFPs.o: firstFPs.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
ForwardSampler.o: ForwardSampler.cpp ForwardSampler.h LPSolver.h \
  Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Random.h
FPListIntersection.o: FPListIntersection.cpp Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h
fplottify.o: fplottify.cpp Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h
generateLinProg.o: generateLinProg.cpp Matrix.h Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h Random.h
generateSparseLinProg.o: generateSparseLinProg.cpp Matrix.h Vector.h \
  Array.h Streamable.h Set.h Array.cpp Numerical.h Random.h
GramSchmidt.o: GramSchmidt.cpp GramSchmidt.h Matrix.h Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h
GridEmissionModel.o: GridEmissionModel.cpp GridEmissionModel.h \
  EmissionModel.h Combinatorics.h
GridLP.o: GridLP.cpp GridLP.h ForwardSampler.h LPSolver.h Numerical.h \
  GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp \
  Random.h
GridSimplex.o: GridSimplex.cpp GridSimplex.h Simplex.h LPSolver.h \
  Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Random.h
GroupPowerBigraph.o: GroupPowerBigraph.cpp GroupPowerBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h BasicGroupBigraph.h \
  ReplicateIndexer.h ReplicateIndexer.cpp BasicBigraph.h Cache.h \
  PivdoSplitter.h PeptideFeatureReader.h KDMethod.h
HashTable.o: HashTable.cpp
interpolateROC.o: interpolateROC.cpp Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h ROCInterpolator.h
KDMethod.o: KDMethod.cpp KDMethod.h Matrix.h Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h Random.h
LPBigraph.o: LPBigraph.cpp LPBigraph.h BasicBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
LPSolver.o: LPSolver.cpp LPSolver.h Numerical.h GramSchmidt.h Matrix.h \
  Vector.h Array.h Streamable.h Set.h Array.cpp Random.h
mainAffineScalingArray.o: mainAffineScalingArray.cpp AffineScaling.h \
  LPSolver.h Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h \
  Streamable.h Set.h Array.cpp Random.h
mainAffineScaling.o: mainAffineScaling.cpp AffineScaling.h LPSolver.h \
  Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Random.h
main.o: main.cpp GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h \
  Set.h Array.cpp Numerical.h
mainCrossValSplit.o: mainCrossValSplit.cpp PivdoSplitter.h BasicBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h
mainEM.o: mainEM.cpp EMBigraph.h BasicGroupBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h ReplicateIndexer.h ReplicateIndexer.cpp \
  BasicBigraph.h Cache.h
mainFido.o: mainFido.cpp GroupPowerBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h BasicGroupBigraph.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h PivdoSplitter.h \
  PeptideFeatureReader.h KDMethod.h
mainFidoDetectability.o: mainFidoDetectability.cpp GroupPowerBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h BasicGroupBigraph.h \
  ReplicateIndexer.h ReplicateIndexer.cpp BasicBigraph.h Cache.h \
  PivdoSplitter.h PeptideFeatureReader.h KDMethod.h
mainFido_Likelihood.o: mainFido_Likelihood.cpp GroupPowerBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h BasicGroupBigraph.h \
  ReplicateIndexer.h ReplicateIndexer.cpp BasicBigraph.h Cache.h \
  PivdoSplitter.h PeptideFeatureReader.h KDMethod.h
mainFidoMCMC.o: mainFidoMCMC.cpp GroupPowerBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h BasicGroupBigraph.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h PivdoSplitter.h \
  PeptideFeatureReader.h KDMethod.h
mainForwardSamplerArray.o: mainForwardSamplerArray.cpp ForwardSampler.h \
  LPSolver.h Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h \
  Streamable.h Set.h Array.cpp Random.h
mainForwardSampler.o: mainForwardSampler.cpp ForwardSampler.h LPSolver.h \
  Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Random.h
mainGridLP.o: mainGridLP.cpp GridLP.h ForwardSampler.h LPSolver.h \
  Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Random.h
mainGridSimplex.o: mainGridSimplex.cpp GridSimplex.h Simplex.h LPSolver.h \
  Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Random.h
mainMCMC2Pivdo.o: mainMCMC2Pivdo.cpp PivdoSplitter.h BasicBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h
mainPivdo.o: mainPivdo.cpp PivdoBigraph.h BasicGroupBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h
mainRandomSimplex.o: mainRandomSimplex.cpp Simplex.h LPSolver.h \
  Numerical.h GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Random.h
mainROC.o: mainROC.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
mainROC_ProtProph.o: mainROC_ProtProph.cpp Vector.h Array.h Streamable.h \
  Set.h Array.cpp Numerical.h
mainSimplex.o: mainSimplex.cpp Simplex.h LPSolver.h Numerical.h \
  GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp \
  Random.h
Matrix.o: Matrix.cpp Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h
meanStdDev.o: meanStdDev.cpp Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h
MSE.o: MSE.cpp Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
Numerical.o: Numerical.cpp Numerical.h
PeptideFeatureReader.o: PeptideFeatureReader.cpp PeptideFeatureReader.h \
  Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp Numerical.h \
  Random.h KDMethod.h
PivdoBigraph.o: PivdoBigraph.cpp PivdoBigraph.h BasicGroupBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h
PivdoSplitter.o: PivdoSplitter.cpp PivdoSplitter.h BasicBigraph.h \
  ProteinIdentifier.h StringTable.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h
plottify.o: plottify.cpp Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h
ProteinBigraph.o: ProteinBigraph.cpp ProteinBigraph.h StringTable.h \
  HashTable.h Array.h Streamable.h Set.h Array.cpp HashTable.cpp \
  StringFunctions.h Matrix.h Vector.h Numerical.h Random.h
ProteinIdentifier.o: ProteinIdentifier.cpp ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
Proteomics.o: Proteomics.cpp Proteomics.h Array.h Streamable.h Set.h \
  Array.cpp BasicBigraph.h ProteinIdentifier.h StringTable.h HashTable.h \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
Random.o: Random.cpp Random.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
rank.o: rank.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
relativeOptimize.o: relativeOptimize.cpp Vector.h Array.h Streamable.h \
  Set.h Array.cpp Numerical.h
ReplicateIndexer.o: ReplicateIndexer.cpp
rlistFDR.o: rlistFDR.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
rlistFPLists.o: rlistFPLists.cpp Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h
rlistROC.o: rlistROC.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
rlistSensitivityAtFDR.o: rlistSensitivityAtFDR.cpp Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h
roc50.o: roc50.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
rocDifferenceInterpolate.o: rocDifferenceInterpolate.cpp Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h ROCInterpolator.h
ROCInterpolator.o: ROCInterpolator.cpp ROCInterpolator.h Array.h \
  Streamable.h Set.h Array.cpp
Set.o: Set.cpp Array.h Streamable.h Set.h Array.cpp
Simplex.o: Simplex.cpp Simplex.h LPSolver.h Numerical.h GramSchmidt.h \
  Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp Random.h
Streamable.o: Streamable.cpp Streamable.h
stringCounter.o: stringCounter.cpp Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h StringTable.h HashTable.h HashTable.cpp \
  StringFunctions.h
testFasta.o: testFasta.cpp FastaReader.h Array.h Streamable.h Set.h \
  Array.cpp Proteomics.h BasicBigraph.h ProteinIdentifier.h StringTable.h \
  HashTable.h HashTable.cpp StringFunctions.h Matrix.h Vector.h \
  Numerical.h Random.h Model.h Combinatorics.h
tmp.o: tmp.cpp Cache.h Model.h Combinatorics.h
Vector.o: Vector.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
w_test.o: w_test.cpp Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
AffineScaling.o: AffineScaling.h LPSolver.h Numerical.h GramSchmidt.h \
  Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp Random.h
Array.o: Array.h Streamable.h Set.h Array.h Array.cpp
BasicBigraph.o: BasicBigraph.h ProteinIdentifier.h StringTable.h \
  HashTable.h Array.h Streamable.h Set.h Array.cpp HashTable.cpp \
  StringFunctions.h Matrix.h Vector.h Numerical.h Random.h Model.h \
  Combinatorics.h
BasicGroupBigraph.o: BasicGroupBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h ReplicateIndexer.h ReplicateIndexer.cpp \
  BasicBigraph.h Cache.h
Cache.o: Cache.h
Combinatorics.o: Combinatorics.h
EMBigraph.o: EMBigraph.h BasicGroupBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h ReplicateIndexer.h ReplicateIndexer.cpp \
  BasicBigraph.h Cache.h
EmissionModel.o: EmissionModel.h Combinatorics.h
FastaReader.o: FastaReader.h Array.h Streamable.h Set.h Array.cpp
ForwardSampler.o: ForwardSampler.h LPSolver.h Numerical.h GramSchmidt.h \
  Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp Random.h
GramSchmidt.o: GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h
GridEmissionModel.o: GridEmissionModel.h EmissionModel.h Combinatorics.h
GridLP.o: GridLP.h ForwardSampler.h LPSolver.h Numerical.h GramSchmidt.h \
  Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp Random.h
GridSimplex.o: GridSimplex.h Simplex.h LPSolver.h Numerical.h \
  GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp \
  Random.h
GroupPowerBigraph.o: GroupPowerBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h BasicGroupBigraph.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h PivdoSplitter.h \
  PeptideFeatureReader.h KDMethod.h
HashTable.o: HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp
KDMethod.o: KDMethod.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h Random.h
LPBigraph.o: LPBigraph.h BasicBigraph.h ProteinIdentifier.h StringTable.h \
  HashTable.h Array.h Streamable.h Set.h Array.cpp HashTable.cpp \
  StringFunctions.h Matrix.h Vector.h Numerical.h Random.h Model.h \
  Combinatorics.h
LPSolver.o: LPSolver.h Numerical.h GramSchmidt.h Matrix.h Vector.h \
  Array.h Streamable.h Set.h Array.cpp Random.h
Matrix.o: Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
Model.o: Model.h Combinatorics.h
Numerical.o: Numerical.h
PeptideFeatureReader.o: PeptideFeatureReader.h Matrix.h Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h Random.h KDMethod.h
PivdoBigraph.o: PivdoBigraph.h BasicGroupBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h ReplicateIndexer.h ReplicateIndexer.cpp \
  BasicBigraph.h Cache.h
PivdoSplitter.o: PivdoSplitter.h BasicBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
ProteinBigraph.o: ProteinBigraph.h StringTable.h HashTable.h Array.h \
  Streamable.h Set.h Array.cpp HashTable.cpp StringFunctions.h Matrix.h \
  Vector.h Numerical.h Random.h
ProteinIdentifier.o: ProteinIdentifier.h StringTable.h HashTable.h \
  Array.h Streamable.h Set.h Array.cpp HashTable.cpp StringFunctions.h \
  Matrix.h Vector.h Numerical.h Random.h Model.h Combinatorics.h
Proteomics.o: Proteomics.h Array.h Streamable.h Set.h Array.cpp \
  BasicBigraph.h ProteinIdentifier.h StringTable.h HashTable.h \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
Random.o: Random.h Array.h Streamable.h Set.h Array.cpp Numerical.h
ReplicateIndexer.o: ReplicateIndexer.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp ReplicateIndexer.cpp
ROCInterpolator.o: ROCInterpolator.h Array.h Streamable.h Set.h Array.cpp
Set.o: Set.h Array.h Streamable.h Set.h Array.cpp
Simplex.o: Simplex.h LPSolver.h Numerical.h GramSchmidt.h Matrix.h \
  Vector.h Array.h Streamable.h Set.h Array.cpp Random.h
Streamable.o: Streamable.h
StringFunctions.o: StringFunctions.h
StringTable.o: StringTable.h HashTable.h Array.h Streamable.h Set.h \
  Array.cpp HashTable.cpp StringFunctions.h
Vector.o: Vector.h Array.h Streamable.h Set.h Array.cpp Numerical.h


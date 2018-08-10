CPPFLAGS=-Wall -ffast-math -pipe -O4 -g
BINPATH=.
VPATH=../../bin
TARGETS=Fido FidoChooseParameters

all: $(TARGETS)

depend:
	g++ -MM *.cpp *.h

clean:
	rm *.o;


Fido: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFido.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

FidoChooseParameters: Set.o Vector.o Numerical.o Matrix.o GramSchmidt.o LPSolver.o ForwardSampler.o Random.o ProteinIdentifier.o BasicBigraph.o BasicGroupBigraph.o GroupPowerBigraph.o PeptideFeatureReader.o KDMethod.o mainFido_ChooseParameters.o PivdoSplitter.o
	g++ $(CPPFLAGS) $^ -o $(BINPATH)/$@

%.o: %.cpp
	g++ $(CPPFLAGS) -c $<

#automade
Array.o: Array.h Streamable.h Set.h Array.h Array.cpp
BasicBigraph.o: BasicBigraph.cpp BasicBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
BasicGroupBigraph.o: BasicGroupBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h ReplicateIndexer.h ReplicateIndexer.cpp \
  BasicBigraph.h Cache.h
Cache.o: Cache.h
Combinatorics.o: Combinatorics.h
ForwardSampler.o: ForwardSampler.h LPSolver.h Numerical.h GramSchmidt.h \
  Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp Random.h
GramSchmidt.o: GramSchmidt.h Matrix.h Vector.h Array.h Streamable.h Set.h \
  Array.cpp Numerical.h
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
LPSolver.o: LPSolver.h Numerical.h GramSchmidt.h Matrix.h Vector.h \
  Array.h Streamable.h Set.h Array.cpp Random.h
mainFido.o: mainFido.cpp GroupPowerBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h BasicGroupBigraph.h ReplicateIndexer.h \
  ReplicateIndexer.cpp BasicBigraph.h Cache.h PivdoSplitter.h \
  PeptideFeatureReader.h KDMethod.h
Matrix.o: Matrix.h Vector.h Array.h Streamable.h Set.h Array.cpp \
  Numerical.h
Model.o: Model.h Combinatorics.h
Numerical.o: Numerical.cpp Numerical.h
PeptideFeatureReader.o: PeptideFeatureReader.h Matrix.h Vector.h Array.h \
  Streamable.h Set.h Array.cpp Numerical.h Random.h KDMethod.h
PivdoSplitter.o: PivdoSplitter.h BasicBigraph.h ProteinIdentifier.h \
  StringTable.h HashTable.h Array.h Streamable.h Set.h Array.cpp \
  HashTable.cpp StringFunctions.h Matrix.h Vector.h Numerical.h Random.h \
  Model.h Combinatorics.h
ProteinIdentifier.o: ProteinIdentifier.h StringTable.h HashTable.h \
  Array.h Streamable.h Set.h Array.cpp HashTable.cpp StringFunctions.h \
  Matrix.h Vector.h Numerical.h Random.h Model.h Combinatorics.h
Random.o: Random.h Array.h Streamable.h Set.h Array.cpp Numerical.h
ReplicateIndexer.o: ReplicateIndexer.h HashTable.h Array.h Streamable.h \
  Set.h Array.cpp HashTable.cpp ReplicateIndexer.cpp
Set.o: Set.h Array.h Streamable.h Set.h Array.cpp
Streamable.o: Streamable.h
StringFunctions.o: StringFunctions.h
StringTable.o: StringTable.h HashTable.h Array.h Streamable.h Set.h \
  Array.cpp HashTable.cpp StringFunctions.h
Vector.o: Vector.h Array.h Streamable.h Set.h Array.cpp Numerical.h


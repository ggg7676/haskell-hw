module Group2_fn45509 where
    
data Base = C | G | A | T
    deriving (Eq, Show)

data GeneticCode = GeneticCode [Base]
    deriving (Eq, Show)

data Codon = Codon Base Base Base
    deriving (Eq, Show)

data Aminoacid = Aminoacid Int
    deriving (Eq, Ord, Show)

type CodingTable = [(Codon, Aminoacid)]
type Gene = [Codon]
type Protein = [Aminoacid]


allCodons :: [Codon]
allCodons =  [
          Codon T T T, Codon T C T, Codon T A T, Codon T G T, Codon T T C, Codon T C C, Codon T A C, Codon T G C, Codon T T A, Codon T C A, Codon T A A, 
          Codon T G A, Codon T T G, Codon T C G, Codon T A G, Codon T G G, Codon C T T, Codon C C T, Codon C A T, Codon C G T, Codon C T C, Codon C C C, 
          Codon C A C, Codon C G C, Codon C T A, Codon C C A, Codon C A A, Codon C G A, Codon C T G, Codon C C G, Codon C A G, Codon C G G, Codon A T T,
          Codon A C T, Codon A A T, Codon A G T, Codon A T C, Codon A C C, Codon A A C, Codon A G C, Codon A T A, Codon A C A, Codon A A A, Codon A G A,
          Codon A T G, Codon A C G, Codon A A G, Codon A G G, Codon G T T, Codon G C T, Codon G A T, Codon G G T, Codon G T C, Codon G C C, Codon G A C, 
          Codon G G C, Codon G T A, Codon G C A, Codon G A A, Codon G G A, Codon G T G, Codon G C G, Codon G A G, Codon G G G
             ]


getCodons :: GeneticCode -> [Codon]
getCodons (GeneticCode []) = []
getCodons (GeneticCode (x:y:z:rest)) = (Codon x y z) : (getCodons (GeneticCode rest))

isStopCodon :: CodingTable -> Codon -> Bool
isStopCodon [] _ = True
isStopCodon ((c, a) : rest) codon = if codon == c then False else isStopCodon rest codon

getGenes :: CodingTable -> [Codon] -> [Gene] 
getGenes _ [] = []
getGenes table list@(codon : next)  = let newList = helper table list
                                          len = (length newList) + 1
                                      in newList : getGenes table (drop len list) where 
                                          helper :: CodingTable -> [Codon] -> Gene
                                          helper _ [] = []
                                          helper table (codon:rest) = if isStopCodon table codon then [] else codon : helper table rest 

getAminoacid :: CodingTable -> Codon -> Aminoacid
getAminoacid ((c,a) : rest) codon = if codon == c then a else getAminoacid rest codon

getProteins :: CodingTable -> [Gene] -> [Protein]
getProteins _ [] = []
getProteins table (codonList:rest) = (helperList table codonList) : (getProteins table rest) where
    helperList :: CodingTable -> [Codon] -> [Aminoacid]
    helperList _ [] = [] 
    helperList table (codon : next) = (getAminoacid table codon) : (helperList table next) 

hasRepeatingElements :: [Protein] -> Bool
hasRepeatingElements [] = False 
hasRepeatingElements (p:rest) = if p `elem` rest then True else hasRepeatingElements rest 


hasSameProteins :: CodingTable -> GeneticCode -> Bool
hasSameProteins table genCode = hasRepeatingElements (getProteins table (getGenes table (getCodons genCode))) 



getCodonsForAminoacid :: CodingTable -> Aminoacid -> [Codon]
getCodonsForAminoacid table aminoacid = [c | (c,a) <- table, a == aminoacid]

getDifferencesList :: Codon -> [Codon] -> [Int]
getDifferencesList _ [] = []
getDifferencesList codon (c:rest) = countDifferences codon c : getDifferencesList codon rest where 
    countDifferences :: Codon -> Codon -> Int 
    countDifferences c1 c2 = let baseList1 = (\(Codon x y z) -> x:y:z:[]) c1
                                 baseList2 = (\(Codon x y z) -> x:y:z:[]) c2
                             in  helper baseList1 baseList2 where 
                                 helper [] [] = 0
                                 helper (b1:rest1) (b2:rest2) = (if b1 == b2 then 0 else 1) + helper rest1 rest2 
                            
                             
countMax :: CodingTable -> [Codon] -> Int 
countMax _ [] = 0
countMax table (currCodon : next) =  let codonsForAminoacid = if isStopCodon table currCodon then [c | c <- allCodons, isStopCodon table c]
                                                                                          else getCodonsForAminoacid table $ getAminoacid table currCodon
                                         diffList = getDifferencesList currCodon codonsForAminoacid
                                     in maximum diffList + countMax table next 

                                                
maxMutations :: CodingTable -> GeneticCode -> Int
maxMutations table g = countMax table $ getCodons g

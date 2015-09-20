public class ProteinTranslator2 {
	
	
	public static void main (String[] args) {
		
		String test = "GCGCGCAGCCGGGATTCCCGCTTGCACGGTTGAAAATGGTTGTCGAACAAGAATTCGCTCAGATCAAACATGTTCTGCATGGTATCAGCCTGCTGGGTCAGTGCCCGGATAGCATCAACGCCGCGCTGATTTGCCGTGGCGAAAAAATGTCGATCGCGATTATGGCGGGACTTCTGGAGGCGCGTGGGCATCGCGTCACGGTGATCGATCCGGTAGAAAAATTGCTGGCGGTGGGCCATTACCTTGAATCTACCGTTGATATCGCGGAATCGACTCGCCGTATCGCCGCCAGCCAGATCCCGGCCGATCACATGATCCTGATGGCGGGCTTTACCGCCGGTAATGAAAAGGGTGAACTGGTGGTGCTGGGCCGTAATGGTTCCGACTATTCCGCCGCCGTGCTGGCCGCCTGTTTACGCGCTGACTGCTGTGAAATCTGGACTGACGTCGATGGCGTGTATACCTGTGACCCGCGCCAGGTGCCGGACGCCAGGCTGCTGAAATCGATGTCCTACCAGGAAGCGATGGAACTCTCTTACTTCGGCGCCAAAGTCCTTCACCCTCGCACCATTACGCCCATCGCCCAGTTCCAGATCCCCTGTCTGATTAAAAATACCGGTAATCCGCAGGCGCCAGGAACGCTGATCGGCGCGTCCAGCGACGATGATAACCTGCCGGTTAAAGGGATCTCTAACCTTAACAACATGGCGATGTTTAGCGTCTCCGGCCCGGGAATGAAAGGGATGATTGGGATGGCGGCGCGTGTTTTCGCCGCCATGTCTCGCGCCGGGATCTCGGTGGTGCTCATTACCCAGTCCTCCTCTGAGTACAGCATCAGTTTCTGTGTGCCGCAGAGTGACTGCGCGCG";
		for (String i : findAllProteins(test)) System.out.println(i);
		
	} //end main method
	
	
	public static String findAminoAcid (String s){
        
        final String[] AA_CODE ={"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"};
        final String[] AA_CODE_VALUE = {"F","F","L","L","S","S","S","S","Y","Y","STOP","STOP","C","C","STOP","W","L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"};
        
        for (int i=0; i < AA_CODE.length; i++){
              if (AA_CODE[i].equals(s)){
                    System.out.println(AA_CODE[i] + "   " + AA_CODE_VALUE[i]);
                    return AA_CODE_VALUE[i];
              }
        }
        
        return "";  
	} //end findAminoAcid
  
	public static String findRNAStrand (String s){
        
        String result = "";
        
        for(int i=0; i<s.length(); i++){
              if(s.charAt(i) == 'A')
                    result = result + 'A';
              if(s.charAt(i) == 'T')
                    result = result + 'U';
              if(s.charAt(i) == 'G')
                    result = result + 'G';
              if(s.charAt(i) == 'C')
                    result = result + 'C';
        }
        
        return result;
  } //end findRNAStrand

	
  public static String findComplement (String s){
        
        String result = "";
        
        for(int i=0; i<s.length(); i++){
              if(s.charAt(i) == 'A')
                    result = result + 'T';
              if(s.charAt(i) == 'T')
                    result = result + 'A';
              if(s.charAt(i) == 'G')
                    result = result + 'C';
              if(s.charAt(i) == 'C')
                    result = result + 'G';
        }
        
        return result;
  } //end findComplement

	
	/** 
	 *  Takes one half of a DNA double helix as String s and
	 *  returns possible protein that can be constructed from the entire
	 *  strand reading from the first nucleotide
	 *  @param  s represents one half of a DNA double helix
	 *  @return a string representing the protein that could be  
	 *  constructed from String s.  If strand does not define a 
	 *  protein return ìINVALID ñ NO STOP CODONî   
	 *  precondition s is a string consisting only of characters A,C,T and  
	 *  G
	 *  
	 **/
	public static String findProteinA(String s) {
		
		final String[] AA_CODE ={"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"};
		final String[] AA_CODE_VALUE = {"F","F","L","L","S","S","S","S","Y","Y","STOP","STOP","C","C","STOP","W","L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"};
		
		s = trimString(s);
		String returning = "";
		
		for (int i = 0; i <= s.length() - 3; i += 3) {
			for (int j = 0; j < AA_CODE.length; j++) {
				if (s.substring(i, i + 3).equals(AA_CODE[j])) {
					if (AA_CODE_VALUE[j].equals("STOP")) return returning;
					else returning += AA_CODE_VALUE[j];
				}
			}
		}
		
		return "INVALID - NO STOP CODON";
	} //end findProteinA

	
	/**
	* Takes one half of a DNA double helix as String s and returns  
	* an array of Strings of length three.  Each element in the 
	* array of Strings represents the constructed protein, assuming a 
	* valid stop codon is found in the reading frame.  If a 
	* reading frame does not have valid stop codon, the array stores 
	* ìINVALID ñ NO STOP CODONî in the array
	* @param String s representing the DNA sequence
	* @return String[] of length 3
	* 
	*/ 
	public static String[] findProteinB(String s) {
		
		String[] returning = new String[3];
		String set = s;
		s = trimString(set);
		String s2 = trimString(set.substring(1));
		String s3 = trimString(set.substring(2));
		
		returning[0] = findProteinA(s);
		returning[1] = findProteinA(s2);
		returning[2] = findProteinA(s3);
		
		return returning;
		
	} //end findProteinB
	
	
	
	/*
	* Takes one half of a DNA double helix as a String s and returns an 
	* array which stores the index of the first nucleotide for all 
	* start codons. 
	* @param String s representing one half of a DNA double helix
	* @return int[] storing the position of the first nucleotide for all 
	* possible start codons. 
	*
	*/
	public static int[] getStartCodonPositions(String s) {
		ArrayList<Integer> prereturn = new ArrayList<Integer>();
		String set = s;
		set = trimString(set);
		for (int i = 0; i <= s.length() - 3; i ++) {
			if (s.substring(i, i + 3).equals("ATG")) prereturn.add(i);
		}
		
		return convertArray(prereturn);
		
	} //end getStartCodonPositions
	
	
	/*
	* Takes one half of a DNA double helix as a String s and returns an 
	* array which stores the index of the first nucleotide for all 
	* stop codons. 
	* @param String s representing one half of a DNA double helix
	* @return int[] storing the position of the first nucleotide for all 
	* possible stop codons. 
	*
	*/
	public static int[] getStopCodonPositions(String s) {
		ArrayList<Integer> prereturn = new ArrayList<Integer>();
		String set = s;
		set = trimString(set);
		for (int i = 0; i <= s.length() - 3; i ++) {
			if (s.substring(i, i + 3).equals("TAA") || s.substring(i, i + 3).equals("TAG") || s.substring(i, i + 3).equals("TGA")) prereturn.add(i);
		}
		
		return convertArray(prereturn);
	} //end getStopCodonPositions

	
	// This method takes one half of a DNA double helix and returns all 
	// possible proteins that can be constructed using any of the 
	// three valid reading frames and containing a start codon.
	// It should return an array of Strings which holds all valid 
	// proteins.
	public static String[] findAllProteins(String s) {
		
		ArrayList<String> store = new ArrayList<String>();
		String r1 = trimString(s);
		String r2 = trimString(s.substring(1));
		String r3 = trimString(s.substring(2));
		
		int[] start1 = getStart(r1);
		int[] stop1 = getStop(r1);
		
		int[] start2 = getStart(r2);
		int[] stop2 = getStop(r2);
		
		int[] start3 = getStart(r3);
		int[] stop3 = getStop(r3);
		
		if (start1.length != 0 && stop1.length != 0) {
			for (int i = 0; i < stop1.length; i++) {
				for (int u = 0; u < start1.length; u++) {
					if (start1[u] < stop1[i]) {
						store.add(findProteinA(r1.substring(start1[u], stop1[i] + 3)));
						start1[u] = 999999999;
					} //while
				} //u
			} //i
		}
		
		if (start2.length != 0 && stop2.length != 0) {
			for (int i = 0; i < stop2.length; i++) {
				for (int u = 0; u < start2.length; u++) {
					if (start2[u] < stop2[i]) {
						store.add(findProteinA(r2.substring(start2[u], stop2[i] + 3)));
						start2[u] = 999999999;
					} //while
				} //u
			} //i
		}
		
		if (start3.length != 0 && stop3.length != 0) {
			for (int i = 0; i < stop3.length; i++) {
				for (int u = 0; u < start3.length; u++) {
					if (start3[u] < stop3[i]) {
						store.add(findProteinA(r3.substring(start3[u], stop3[i] + 3)));
						start3[u] = 999999999;
					} //while
				} //u
			} //i
		}
		
		return convertStringArray(store);
		
	} //end findAllProteins

	
	
	
	
	//Tool Classes:
	
	/**
	 * Trims string so that it is divisible by 3
	 * Helps to eliminate unnecessary letters
	 * @param s
	 * @return
	 */
	public static String trimString(String s) {
		
		while (s.length() % 3 != 0) {
			s = s.substring(0, s.length() - 1);
		}
		return s;
	} //end trimString
	
	/**
	 * Simple method used to convert an Integer ArrayList to an int array
	 * @param al
	 * @return
	 */
	public static int[] convertArray(ArrayList<Integer> al) {
		int[] returning = new int[al.size()];
		for (int i = 0; i < returning.length; i++) returning[i] = al.get(i);
		return returning;
	} //end convertArray
	
	/**
	 * Convert a String ArrayList to a string array
	 * @param al
	 * @return
	 */
	public static String[] convertStringArray(ArrayList<String> al) {
		String[] returning = new String[al.size()];
		for (int i = 0; i < returning.length; i++) returning[i] = al.get(i);
		return returning;
	} //end convertStringArray
	
	/**
	 * Alternate method to getStartCodons because getStartCodons does not go by reading frames
	 * @param s
	 * @return
	 */
	public static int[] getStart(String s) {
		ArrayList<Integer> prereturn = new ArrayList<Integer>();
		String set = s;
		set = trimString(set);
		for (int i = 0; i <= s.length() - 3; i += 3) {
			if (s.substring(i, i + 3).equals("ATG")) prereturn.add(i);
		}
		
		return convertArray(prereturn);
		
	} //end getStartCodonPositions
	
	/**
	 * Alternate method to getStopCodons, this one counts with reading frames in mind
	 * @param s
	 * @return
	 */
	public static int[] getStop(String s) {
		ArrayList<Integer> prereturn = new ArrayList<Integer>();
		String set = s;
		set = trimString(set);
		for (int i = 0; i <= s.length() - 3; i += 3) {
			if (s.substring(i, i + 3).equals("TAA") || s.substring(i, i + 3).equals("TAG") || s.substring(i, i + 3).equals("TGA")) prereturn.add(i);
		}
		
		return convertArray(prereturn);
	} //end getStopCodonPositions

	
} //end class

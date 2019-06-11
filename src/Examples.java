import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class Examples {

	private final static StringBuilder SEQ1 = new StringBuilder("ATTGCCATT");
	private final static StringBuilder SEQ2 = new StringBuilder("ATGGCCATT");
	private final static StringBuilder SEQ3 = new StringBuilder("ATCCAATTTT");
	private final static StringBuilder SEQ4 = new StringBuilder("ATCTTCTT");
	private final static StringBuilder SEQ5 = new StringBuilder("ACTGACC");

	public static void main(String[] args) {

		final List<StringBuilder> seqs = new ArrayList<StringBuilder>();
		seqs.add(SEQ2);
		seqs.add(SEQ3);
		seqs.add(SEQ4);
		seqs.add(SEQ5);

		(starAlignment(SEQ1, seqs)).forEach(System.out::println);
	}
	
	private static double shannonEntropy(String s) {
		return s.chars().mapToDouble(i -> -Math.log((double) s.codePoints()
				.filter(c -> c == i).count() / s.length())
				/ Math.log(2) / s.length()).sum();
	}

	private static float[][] multipleAlignment(ArrayList<StringBuilder> seqs1, ArrayList<StringBuilder> seqs2) {

		float[][] dp = new float[seqs2.get(0).length() + 1][seqs1.get(0).length() + 1];

		for (int i = 0; i <= seqs2.get(0).length(); i++)
			dp[i][0] = i * Score.MSA_GAP;

		for (int j = 1; j <= seqs1.get(0).length(); j++)
			dp[0][j] = j * Score.MSA_GAP;
		
		List<Character> chars2 = new ArrayList<>();
		List<Character> chars1 = new ArrayList<>();

		for (int i = 1; i <= seqs2.get(0).length(); i++) {
			for (int j = 1; j <= seqs1.get(0).length(); j++) {
				
				chars2 = charsByIndex(i - 1, seqs2);
				chars1 = charsByIndex(j - 1, seqs1);
				chars1.addAll(chars2);
				
				float scoreDiag = dp[i - 1][j - 1]
						+ (sumOfPairs(chars1));
				
				float scoreLeft = dp[i][j - 1] + Score.MSA_GAP;
				float scoreUp = dp[i - 1][j] + Score.MSA_GAP;

				dp[i][j] = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
			}
		}

		return dp;
	}

	private static List<Character> charsByIndex(int index, ArrayList<StringBuilder> seqs) {

		List<Character> chars = new ArrayList<>();

		for (StringBuilder seq : seqs) {
			chars.add(seq.charAt(index));
		}

		return chars;
	}

	private static float sumOfPairs(List<Character> chars) {
		float score = 0.0f;

		Iterator<Character> i = chars.iterator();
		while (i.hasNext()) {
			Character c = i.next();
			
			i.remove();
			
			Iterator<Character> i2 = chars.iterator();
			while (i2.hasNext()) {				
				if (c == i2.next()) {
					score = score + Score.MSA_MATCH;
				} else {
					score = score + Score.MSA_MISMATCH;
				}
			}
			
			
		}

		return score;
	}

	/**
	 * Merging the sequences in star alignment:
	 * • Use the anchor as the "guide" sequence
	 * • Add iteratively each pair-wise alignment to the multiple alignment
	 * • Go column by column:
	 * – If there is no gap neither in the guide sequence in the multiple
	 * alignment nor in the merged alignment (or both have gaps)
	 * simply put the letter paired with the guide sequence into the
	 * appropriate column (all steps of the first merge are of this type.
	 * – If pair-wise alignment produced a gap in the guide sequence,
	 * force the gap on the whole column of already aligned sequences
	 * (compare second merge).
	 * – If there is a gap in added sequence but not in the guide
	 * sequences, keep the gap in the added sequence.
	 * 
	 * @param anchor Anchor
	 * @param seqs Sequences
	 * 
	 * @return result
	 */
	private static List<StringBuilder> starAlignment(final StringBuilder anchor, final List<StringBuilder> seqs) {

		List<StringBuilder> anchors = new ArrayList<StringBuilder>();

		// align
		for (int i = 0; i < seqs.size(); i++) {
			int[][] dp = align(anchor, seqs.get(i));
			StringBuilder[] aligned = tracebackDownmost(dp, anchor, seqs.get(i));
			anchors.add(aligned[0]);
			seqs.set(i, aligned[1]);
		}

		StringBuilder updatedAnchor = new StringBuilder(anchor);

		// updatedAnchor
		for (int i = 0; i < anchors.size(); i++) {
			int[][] dp = align(updatedAnchor, anchors.get(i));
			StringBuilder[] aligned = tracebackDownmost(dp, updatedAnchor, anchors.get(i));
			updatedAnchor = aligned[0];
		}

		List<StringBuilder> result = new ArrayList<StringBuilder>();
		result.add(updatedAnchor);

		// onceGapAlwaysGap
		for (int i = 0; i < seqs.size(); i++) {
			int[][] dp = align(updatedAnchor, seqs.get(i));
			StringBuilder[] aligned = tracebackDownmost(dp, updatedAnchor, seqs.get(i));
			result.add(aligned[1]);
		}

		return result;
	}

	@SuppressWarnings("unused")
	private static StringBuilder[] tracebackUpmost(int[][] dp,
			StringBuilder seq1,
			StringBuilder seq2) {

		StringBuilder seqA = new StringBuilder();
		StringBuilder seqB = new StringBuilder();

		for (int j = seq1.length(), i = seq2.length(); i > 0 || j > 0;) {

			int scoreDiag = dp[i - 1][j - 1]
					+ (seq2.charAt(i - 1) == seq1.charAt(j - 1) ? Score.MATCH : Score.MISMATCH);

			if (i > 0 && dp[i][j] == dp[i - 1][j] + Score.GAP) {
				seqB.append(seq2.charAt(--i));
				seqA.append("-");
			} else if (i > 0 && j > 0 && dp[i][j] == scoreDiag) {
				seqB.append(seq2.charAt(--i));
				seqA.append(seq1.charAt(--j));
			} else if (j > 0 && dp[i][j] == dp[i][j - 1] + Score.GAP) {
				seqA.append(seq1.charAt(--j));
				seqB.append("-");
			}
		}

		return new StringBuilder[] { seqA.reverse(), seqB.reverse() };
	}

	private static StringBuilder[] tracebackDownmost(int[][] dp,
			StringBuilder seq1,
			StringBuilder seq2) {

		StringBuilder seqA = new StringBuilder();
		StringBuilder seqB = new StringBuilder();

		for (int j = seq1.length(), i = seq2.length(); i > 0 || j > 0;) {

			int scoreDiag = dp[i - 1][j - 1]
					+ (seq2.charAt(i - 1) == seq1.charAt(j - 1) ? Score.MATCH : Score.MISMATCH);

			if (j > 0 && dp[i][j] == dp[i][j - 1] + Score.GAP) {
				seqA.append(seq1.charAt(--j));
				seqB.append("-");
			} else if (i > 0 && j > 0 && dp[i][j] == scoreDiag) {
				seqB.append(seq2.charAt(--i));
				seqA.append(seq1.charAt(--j));
			} else if (i > 0 && dp[i][j] == dp[i - 1][j] + Score.GAP) {
				seqB.append(seq2.charAt(--i));
				seqA.append("-");
			}
		}

		return new StringBuilder[] { seqA.reverse(), seqB.reverse() };
	}

	private static int[][] align(StringBuilder seq1, StringBuilder seq2) {
		int[][] dp = new int[seq2.length() + 1][seq1.length() + 1];

		for (int i = 0; i <= seq2.length(); i++)
			dp[i][0] = i * Score.GAP;

		for (int j = 1; j <= seq1.length(); j++)
			dp[0][j] = j * Score.GAP;

		for (int i = 1; i <= seq2.length(); i++) {
			for (int j = 1; j <= seq1.length(); j++) {
				int scoreDiag = dp[i - 1][j - 1]
						+ (seq2.charAt(i - 1) == seq1.charAt(j - 1) ? Score.MATCH : Score.MISMATCH);
				int scoreLeft = dp[i][j - 1] + Score.GAP;
				int scoreUp = dp[i - 1][j] + Score.GAP;

				dp[i][j] = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
			}
		}

		return dp;
	}
	
}

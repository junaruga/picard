package picard.fingerprint;

import java.util.Arrays;
import java.util.Objects;

public class FingerprintingTestUtils {

    public static boolean areHaplotypeProbabilitiesEqual(final HaplotypeProbabilities lhs, final HaplotypeProbabilities rhs){
        if (lhs == rhs) {
            return true;
        }
        if (lhs == null || rhs == null || lhs.getClass() != rhs.getClass()) {
            return false;
        }

        if (!Objects.equals(lhs.getHaplotype(), rhs.getHaplotype())){
            return false;
        }

        return Arrays.equals(lhs.getLikelihoods(),rhs.getLikelihoods());
    }
}

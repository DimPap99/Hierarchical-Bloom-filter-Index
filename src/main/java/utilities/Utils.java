package utilities;

public class Utils {

    //Doing value >>> <level> is probably faster tho
    public static int getIntervalIndex(int domainTotalBits, int level, int value){
//        int bitsToExtract = totalBits - level;
//        //bring one to the position of MSB for the amount of bits we want to extract
//        //everything else remains 0. -1 to flip all the other bits. This will give us as many '111...' as the
//        //total bits we want to extract. Then shift left by level in order to push them to the leftmost side
//        int mask          = ((1 << bitsToExtract) - 1) << level;
//        //Logical and between our value and mask. The mask will grab the MSB bits we want. Then shift right (divide) based
//        //on the level to find at which interval the value exists
//        //Example: Level 4, max level = 5. 2^5 = 32. So for level 4 we got 32 / 2^4 = 32/16 = 2. We got 2 16 intervals
//        //Consider value = 17. We want to find at which interval it should be placed. Since we got 2 intervals, the total
//        //bits we want to extract are 1, as it can represent 2 options.
//        //Our mask then becomes: 10000. 17: 10001. We do: 17 & 10000 = 10000 (16) and then shift as many positions as the
//        //current level, that is 4. 10000 >> 4 = 1. So it belongs in the second interval.
//
//        return (value & mask) >>> level;
        return value >>> domainTotalBits - level;
    }

    public static enum MemPolicy {
        NONE,
        REACTIVE,
        PREDICTIVE
    }
    /**
     * Formats an integer into  "bitinterval<interval_bits><remaining_bits>".
     *
     * @param totalBits total number of bits in the domain (e.g. 5 for 0..31)
     * @param level     number of binary splits (level = 1 ⇒ 2 intervals, etc.)
     * @param value     value to classify
     * @return          formatted string, e.g. "bitinterval10001"
     */
    /**
     * Returns a string whose left-hand part is the interval index (width = level)
     * and whose right-hand part is a run of '#' characters (count = totalBits-level).
     *
     * Example (totalBits = 5):
     *   value = 17 (10001₂), level = 1  -->  "1####"
     *   value = 22 (10110₂), level = 3  -->  "101##"
     */
    /**
     * Returns  "<interval bits><# placeholders>".
     * Example (totalBits = 5):
     *   value = 22 (10110₂), level = 3  ->  "101##"
     */
    public static String intervalWithHashes(int totalBits, int level, int value) {

        int lowerBits = totalBits - level;          // how many '#' chars
        int interval  = (level == 0) ? 0           // avoid >>> with shift = totalBits
                : value >>> lowerBits;

        /* ---- left part: interval bits ---- */
        String intervalBits;
        if (level > 0) {                            // <-- guard fixes %0s problem
            intervalBits = Integer.toBinaryString(interval);
            intervalBits = String.format("%" + level + "s", intervalBits)
                    .replace(' ', '0');
        } else {
            intervalBits = "";                      // level == 0 ⇒ no bits
        }

        /* ---- right part: '#' placeholders ---- */
        String hashes = "#".repeat(lowerBits);      // repeat(0) returns ""

        return intervalBits + hashes;
    }


}

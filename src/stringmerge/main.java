package stringmerge;

import java.util.Set;

public class main {

    public static void main(String[] args){
        //String output = mergeAlternately("ab", "pqrs");
        //System.out.println(output);

        //String output1 = gcdOfStrings("ABABAB", "ABAB");
        //System.out.println(output1);

//        int[] input2 = new int[]{1, 0, 0, 0, 1};
//        boolean output2 = canPlaceFlowers(input2, 1);
//        System.out.println(output2);

//        int[] input3 = new int[]{0, 1, 0};
//        boolean output3 = canPlaceFlowers(input3, 1);
//        System.out.println(output3);

//        String input = "ai";
//        String reverseVowels = reverseVowels(input);
//        System.out.println(reverseVowels);

        
    }

    public static String mergeAlternately(String word1, String word2) {
        StringBuilder sb = new StringBuilder();
        int word1length = word1.length();
        int word2length = word2.length();
        int maxLength;

        if (word1length > word2length) {
            maxLength = word1length;
        } else {
            maxLength = word2length;
        }

        char[] word1CharArray = word1.toCharArray();
        char[] word2CharArray = word2.toCharArray();

        for(int i = 0; i < maxLength * 2; i++) {
            if (i < word1CharArray.length) {
                sb.append(word1CharArray[i]);
            }
            if (i < word2CharArray.length) {
                sb.append(word2CharArray[i]);
            }
        }

        return sb.toString();
    }

    public static String gcdOfStrings(String str1, String str2) {
        int len1 = str1.length();
        int len2 = str2.length();
        int minLength = Math.min(len1, len2);

        for(int i = minLength; i > 0; i--) {
            if (isDivisor(str1.substring(0, i), str1, str2)) {
                return str1.substring(0, i);
            }
        }

        return "";
    }

    public static boolean isDivisor(String input, String str1, String str2) {
        int len1 = str1.length();
        int len2 = str2.length();
        int inputLength = input.length();

        if (len1 % input.length() != 0 ||  len2 % input.length() != 0 ){
            return false;
        }

        int factor1 = len1 / inputLength;
        int factor2 = len2 / inputLength;

        if (input.repeat(factor1).equals(str1) && input.repeat(factor2).equals(str2)) {
            return true;
        }
        return false;
    }

    public static boolean canPlaceFlowers(int[] flowerbed, int n) {

        int currentIterator = 0;

        if (flowerbed.length == 1 && flowerbed[currentIterator] == 0) {
            return true;
        }

        while (currentIterator < flowerbed.length) {


            if (canPlaceInFlowerBed(flowerbed, currentIterator - 1, currentIterator + 1, currentIterator)) {
                flowerbed[currentIterator] = 1;
                n--;
            }
            currentIterator++;
        }

        return n <= 0;
    }

    public static boolean canPlaceInFlowerBed(int[] flowerbed, int left, int right, int current) {
        if (flowerbed[current] == 0) {

            boolean leftPlotEmpty = (current == 0 || flowerbed[left] == 0);
            boolean rightPlotEmpty = (current == flowerbed.length-1 || flowerbed[right] == 0);

            if (leftPlotEmpty && rightPlotEmpty) {
                return true;
            }
        }

        return false;
    }

    public static String reverseVowels(String s) {
        char[] charArray = s.toCharArray();
        char[] vowels = {'a', 'e', 'i', 'o', 'u', 'A', 'E', 'I', 'O', 'U'};

        int i = 0;
        int j = charArray.length - 1;
        boolean swapleft = false;
        boolean swapright = false;
        while(i <= j) {
            if (containsChar(vowels, charArray[i])) {
                swapleft = true;
            } else {
                i++;
            }

            if (containsChar(vowels, charArray[j])){
                swapright = true;
            } else {
                j--;
            }

            if (swapleft && swapright) {
                char replace = charArray[j];
                charArray[j] = charArray[i];
                charArray[i] = replace;
                swapleft = false;
                swapright = false;
                i++;
                j--;
            }
        }

        return String.valueOf(charArray);
    }

    public static boolean containsChar(char[] array, char target) {
        for (char c : array) {
            if (c == target) {
                return true;
            }
        }
        return false;
    }

}


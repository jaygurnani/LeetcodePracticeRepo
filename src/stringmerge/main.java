package stringmerge;

import com.sun.source.tree.Tree;

import javax.swing.tree.TreeNode;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class main {

    public static void main(String[] args){
//        String output = mergeAlternately("ab", "pqrs");
//        System.out.println(output);

//        String output1 = gcdOfStrings("ABABAB", "ABAB");
//        System.out.println(output1);

//        int[] input2 = new int[]{1, 0, 0, 0, 1};
//        boolean output2 = canPlaceFlowers(input2, 1);
//        System.out.println(output2);

//        int[] input3 = new int[]{0, 1, 0};
//        boolean output3 = canPlaceFlowers(input3, 1);
//        System.out.println(output3);

//        String input = "ai";
//        String reverseVowels = reverseVowels(input);
//        System.out.println(reverseVowels);

//        String input = "a good   example";
//        String output = reverseWords(input);
//        System.out.println(output);

        //int[] input = {0,1,0,3,12};
//        int[] input = {1, 0};
//        moveZeroes(input);

//        int[] input = {3,1,3,4,3};
//        int response = maxOperations2(input, 6);
//        System.out.println(response);

//        String t = "ahbgdc";
//        String s = "abc";
//        boolean output = isSubsequence(s, t);
//        System.out.println(output);

//        int[] nums = {-1};
//        int k = 1;
//        double returnVal = findMaxAverage(nums, k);
//        System.out.println(returnVal);

//        int n = 25;
//        int output = tribonacci(n);
//        System.out.println(output);

//        int[] input = {1,100,1,1,1,100,1,1,100,1};
//        int output = minCostClimbingStairs(input);
//        System.out.println(output);

//        int[] input = {2,7,9,3,1};
//        int output = rob(input);
        int[] input = {3,2,1,5,6,4};
        int k = 2;
        int output = findKthLargest(input, k);
        System.out.println(output);
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

    public static String reverseWords(String s) {
        String after = s.trim().replaceAll(" +", " ");
        char[] charArray = after.toCharArray();
        StringBuilder sb = new StringBuilder();

        int j = 0;
        for(int i = charArray.length - 1; i > -2; i--) {

            if (i == -1 || charArray[i] == ' '){
                String substring = getSubstringFromIndex(charArray, i +1, j);
                sb.append(substring);
                sb.append(" ");
                j = 0;
            } else {
                j++;
            }
        }

       return sb.toString().trim();
    }

    public static String getSubstringFromIndex(char[] charArray, int start, int end) {
        StringBuilder sb = new StringBuilder();
        for(int i = start; i < start + end; i++) {
            sb.append(charArray[i]);
        }

        return sb.toString();
    }

    public static void moveZeroes(int[] nums) {
        int numLength = nums.length;
        int zeroIndex = 0;
        int numberIndex = 0;

        while (zeroIndex != numLength && numberIndex != numLength) {
            boolean containsZero = false;
            boolean containsNumber = false;
            if (nums[zeroIndex] == 0) {
                containsZero = true;
            } else {
                zeroIndex++;
            }

            if (nums[numberIndex] != 0) {
                containsNumber = true;
            } else {
                numberIndex++;
            }

            if (containsZero && containsNumber) {
                if (numberIndex > zeroIndex) {
                    int temp = nums[zeroIndex];
                    nums[zeroIndex] = nums[numberIndex];
                    nums[numberIndex] = temp;
                }

                if (zeroIndex > numberIndex) {
                    numberIndex++;
                }
            }


        }

        System.out.println(Arrays.toString(nums));
    }

    public static int maxOperations(int[] nums, int k) {
        int numberOfOperations = 0;
        ArrayList<Integer> numsList = IntStream.of(nums).boxed().collect(Collectors.toCollection(ArrayList::new));

        int i = 0;
        while (i < numsList.size()) {
            int currentSize = numsList.size();
            numsList = findAndRemoveMatchingK(numsList, i, k);
            if (numsList.size() != currentSize){
                numberOfOperations++;
                i = 0;
            } else {
                i++;
            }

        }

        return numberOfOperations;
    }

    public static ArrayList<Integer> findAndRemoveMatchingK(ArrayList<Integer> nums, int startingIndex, int k) {
        int startingValue = nums.get(startingIndex);

        for (int i = startingIndex+1; i < nums.size() ; i++) {
            int endingValue = nums.get(i);
            if (startingValue + endingValue == k) {
                nums.remove(startingIndex);
                nums.remove(i-1);
                return nums;
            }
        }

        return nums;
    }

    public static int maxOperations2(int nums[], int k) {
        Arrays.sort(nums);

        int left= 0;
        int right = nums.length-1;
        int count= 0;

        while(left < right) {
            int sum = nums[left] + nums[right];

            if (sum == k) {
                left++;
                right--;
                count++;
            } else if (sum > k) {

                right--;
            } else {
                left++;
            }
        }

        return count;
    }

    public static boolean isSubsequence(String s, String t) {
        char[] testString = t.toCharArray();
        int testLength = testString.length;
        int testPointer = 0;

        char[] subString  = s.toCharArray();
        int subLength = subString.length;
        int subPointer = 0;
        boolean result = false;

        if (subLength == 0) {
            return true;
        }

        while (testPointer < testLength && subPointer < subLength) {
            if (subString[subPointer] == testString[testPointer]) {
                subPointer++;
                testPointer++;
            } else {
                testPointer++;
            }
        }

        if (subPointer == subLength && subString[subPointer-1] == testString[testPointer-1]) {
            result = true;
        }

        return result;
    }

    public static int maxDepth(TreeNode root) {
        int left = 0;
        int right = 0;

        if (root == null) {
            return 0;
        }

        if (root.left == null && root.right == null) {
            return 1;
        }

        left = 1 + maxDepth(root.left);
        right =  1 + maxDepth(root.right);

        if (left > right) {
            return left;
        } else {
            return right;
        }
    }

    public boolean leafSimilar(TreeNode root1, TreeNode root2) {
        String root1String = getLeafNode(root1);
        String root2String = getLeafNode(root2);

        return root1String.equals(root2String);
    }

    public String getLeafNode(TreeNode root) {
        String left;
        String right;

        if (root == null) {
            return "";
        }

        if (root.left == null && root.right == null) {
            return String.valueOf(root.val);
        }

        left = getLeafNode(root.left);
        right = getLeafNode(root.right);

        if (left == "") {
            return right;
        }
        if (right == "") {
            return left;
        }

        return left + "," + right;
    }


    public static double findMaxAverage(int[] nums, int k) {
        double currentAverage = Double.NEGATIVE_INFINITY;
        for(int i = 0; i <= nums.length - k; i++) {
            double newAverage = calculateAverage(i, i+k, nums, k);
            if (newAverage > currentAverage) {
                currentAverage = newAverage;
            }
        }

        return currentAverage;
    }

    public static double calculateAverage(int start, int end, int[] nums, int size) {

        int runningTotal = 0;
        for(int i = start; i < end; i++) {
            runningTotal = runningTotal + nums[i];
        }

        return runningTotal / (double) size;
    }

    public static int tribonacci(int n) {
        HashMap<Integer, Integer> dictionary = new HashMap<>();
        dictionary.put(0, 0);
        dictionary.put(1, 1);
        dictionary.put(2, 1);

        for(int i = 0; i <= n; i++) {
            if (dictionary.containsKey(i)) {
                continue;
            } else {
                int toAdd = dictionary.get(i-3) + dictionary.get(i-2) + dictionary.get(i-1);
                dictionary.put(i, toAdd);
            }
        }

        return dictionary.get(n);
    }

    public static int minCostClimbingStairs(int[] cost) {
        HashMap<Integer, Integer> dictionary = new HashMap<>();
        dictionary.put(0, cost[0]);
        dictionary.put(1, cost[1]);

        for (int i = 0; i <= cost.length; i++){
            if (dictionary.containsKey(i)) {
                continue;
            } else {
                int internalCost = 0;
                if (i != cost.length) {
                    internalCost = cost[i];
                }
                int twoStep = dictionary.get(i - 2) + internalCost;
                int oneStep = dictionary.get(i - 1) + internalCost;

                if (twoStep < oneStep) {
                    dictionary.put(i, twoStep);
                } else {
                    dictionary.put(i, oneStep);
                }
            }
        }

        return dictionary.get(cost.length);
    }

    public static int rob(int[] nums) {
        HashMap<Integer, Integer> dictionary = new HashMap<>();
        dictionary.put(0, 0);
        dictionary.put(1, nums[0]);

        for(int i = 2; i < nums.length + 1; i++){
            int singleHouse = dictionary.get(i-1);
            int twoHouses = nums[i-1] + dictionary.get(i-2);
            int max = Math.max(singleHouse, twoHouses);

            dictionary.put(i, max);
        }

        System.out.println(dictionary);
        return dictionary.get(nums.length);
    }

    public static int findKthLargest(int[] nums, int k) {
        PriorityQueue<Integer> pq = new PriorityQueue<>(Collections.reverseOrder());
        int result = 0;

        for(int i = 0; i < nums.length; i++) {
            pq.add(nums[i]);
        }

        for(int i = 0; i < k; i++){
            result = pq.poll();
        }

        return result;
    }

    /* Defined Classes */
    public class TreeNode {
        int val;
        TreeNode left;
        TreeNode right;
        TreeNode() {}
        TreeNode(int val) { this.val = val; }
        TreeNode(int val, TreeNode left, TreeNode right) {
            this.val = val;
            this.left = left;
            this.right = right;
        }
    }
}



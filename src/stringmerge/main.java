package stringmerge;


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
//        int[] input = {3,2,1,5,6,4};
//        int k = 2;
//        int output = findKthLargest(input, k);

        //int output = uniquePaths(3, 7);
        //int output = longestCommonSubsequence("abcde", "ace");

//        int[] input = {1,3,2,8,4,9};
//        int output = maxProfit(input, 3);

//        String s = "abab";
//        String p = "ab";
//        List<Integer> output = findAnagrams(s,p);

        Trie trie = new Trie();
        trie.insert("apple");
        boolean output = trie.search("apple");

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

    public static int uniquePaths(int y, int x) {

        int[][] boardPaths = new int[y][x];
        for (int i = 0; i < y; i++){
            boardPaths[i][0] = 1;
        }
        for (int j = 0; j < x; j++){
            boardPaths[0][j] = 1;
        }


        for(int i = 1; i < y; i++){
            for (int j = 1; j < x; j++) {
                boardPaths[i][j] = boardPaths[i-1][j] + boardPaths[i][j-1];
            }
        }

        return boardPaths[y-1][x-1];

    }

    public static int longestCommonSubsequence(String text1, String text2) {
        int text1Length = text1.toCharArray().length;
        char[] text1CharArray = text1.toCharArray();

        int text2Length = text2.toCharArray().length;
        char[] text2CharArray = text2.toCharArray();

        int[][] dp = new int[text1Length + 1][text2Length + 1];

        for(int i = 1; i < text1Length + 1; i++) {
            for(int j = 1; j < text2Length + 1; j++) {
                if (text1CharArray[i-1] == text2CharArray[j-1]) {
                    dp[i][j] = dp[i-1][j-1] + 1;

                } else {
                    dp[i][j] = Math.max(dp[i-1][j], dp[i][j-1] );
                }
            }
        }

        return dp[text1Length][text2Length];
    }

    public static int maxProfit(int[] prices, int fee) {
        int[][] dp = new int[prices.length][2];

        // if we have no stock. Now we have no stock
        dp[0][0] = 0;
        // if we buy the stock. Now we have the stock
        dp[0][1] = -prices[0];

        for(int i = 1; i < prices.length; i++) {
            // keep the stock or sell it
            dp[i][0] = Math.max(dp[i-1][0], dp[i-1][1] + prices[i] - fee);

            // buy the stock or keep
            dp[i][1] = Math.max(dp[i-1][0] - prices[i], dp[i-1][1]);
        }


        return dp[prices.length-1][0];
    }

    public static ListNode reverseList(ListNode head) {
        if (head == null){
            return null;
        }

        if (head.next == null) {
            return head;
        }

        ListNode curr = head;
        ListNode nextItem = head.next;
        curr.next = null;

        while(nextItem.next != null) {
            ListNode thirdNode = nextItem.next;

            System.out.println("Moving:" + nextItem.next.val + "to" + curr.val);
            nextItem.next = curr;
            curr = nextItem;
            nextItem = thirdNode;
        }
        nextItem.next = curr;

        ListNode print = nextItem;
        while (print.next != null) {
            System.out.println(print.val);
            print = print.next;
        }
        return nextItem;
    }

    // First attempted but timed out
    public static List<Integer> findAnagrams(String s, String p) {
        char[] sChar = s.toCharArray();
        char[] pChar = p.toCharArray();
        List<Integer> result = new ArrayList<>();

        for(int i = 0; i < s.length(); i++) {
            // if we have any character that is contained in the look up
            if (p.contains(String.valueOf(sChar[i]))) {

                // get the first p length characters
                int offset = p.length() + i;
                if (p.length() + i > s.length()) {
                    offset = s.length();
                }
                String substring = s.substring(i, offset);

                if (isAnagram(substring, p)){
                    result.add(i);
                }

            }
        }

        return result;
    }

    public static boolean isAnagram(String s, String p) {
        char[] sChar = s.toCharArray();
        char[] pChar = p.toCharArray();

        for(int i = 0; i < s.length(); i++) {
            for(int j = 0; j < p.length(); j++) {
                if (sChar[i] == pChar[j]) {
                    sChar[i] = 0;
                    pChar[j] = 0;
                }
            }
        }

        boolean sCharIsEmpty = true;
        for(int i = 0; i < sChar.length; i++) {
            if (sChar[i] != 0){
                sCharIsEmpty = false;
                continue;
            }
        }
        boolean pCharIsEmpty = true;
        for(int i = 0; i < pChar.length; i++) {
            if (pChar[i] != 0){
                pCharIsEmpty = false;
                continue;
            }
        }

        if (sCharIsEmpty && pCharIsEmpty) {
            return true;
        }
        return false;
    }

    // Second attempt looking at solution
    public List<Integer> findAnagrams2(String s, String p) {
        int lenS, lenP;
        lenS = s.length();
        lenP = p.length();

        List<Integer> res = new ArrayList<>();
        if(lenP > lenS) return res;

        int[] charsInP = new int[26];
        calculateFreq(p, charsInP);

        int[] charsInS = new int[26];

        // Check the first X character to see if they match
        calculateFreq(s.substring(0, lenP), charsInS);
        if(equalFreq(charsInP, charsInS)){
            res.add(0);
        }

        // Sliding window
        for(int i = lenP; i < lenS; i++){
            char charBeforeStartOfWindow = s.charAt(i-lenP);
            char charAtEndOfWindow = s.charAt(i);

            // Move the Window
            charsInS[charBeforeStartOfWindow - 'a']--;
            charsInS[charAtEndOfWindow - 'a']++;

            // If the frequency count is the same
            if(equalFreq(charsInP, charsInS)) {
                res.add(i - lenP+1);
            }
        }
        return res;
    }

    private void calculateFreq(String s, int[] chars){
        for(char ch : s.toCharArray()){
            chars[ch-'a']++;
        }
    }
    private boolean equalFreq(int[] chars1, int[] chars2){
        for(int i = 0 ; i < 26; i++){
            if (chars1[i] != chars2[i]){
                return false;
            }
        }
        return true;
    }

    public static boolean isValid(String s) {
        char[] sCharArray = s.toCharArray();
        Stack<Character> stack = new Stack<>();

        for(int i = 0; i < s.length(); i++) {
            char charAti = sCharArray[i];

            // We have an input character
            if (charAti == '(' || charAti == '[' || charAti == '{' ){
                stack.push(charAti);
            } else {
                if (stack.size() == 0) {
                    return false;
                }

                // We have an exit character
                char result = stack.pop();
                if (result == '(' && charAti != ')') {
                    return false;
                }
                if (result == '{' && charAti != '}') {
                    return false;
                }
                if (result == '[' && charAti != ']') {
                    return false;
                }
            }
        }

        if (stack.size() == 0) {
            return true;
        }

        return false;
    }

    // Backtracking
    public static List<String> letterCombinations(String digits) {
        if (digits.isEmpty()) return Collections.emptyList();
        String[] phone_map = {"abc", "def", "ghi", "jkl", "mno", "pqrs", "tuv", "wxyz"};
        List<String> output = new ArrayList<>();
        backtrack("", digits, phone_map, output);

        return output;
    }

    public static void backtrack(String combination, String next_digits, String[] phone_map, List<String> output) {
        if (next_digits.isEmpty()) {
            output.add(combination);
        } else {
            int toCheck = next_digits.charAt(0);
            // 0 & 1 don't have any phone items
            String lookupItem = phone_map[toCheck - '2'];
            char[] lookupItemCharArray = lookupItem.toCharArray();

            for (int i = 0; i < lookupItem.length(); i++) {
                backtrack(combination + lookupItemCharArray[i], next_digits.substring(1), phone_map, output );
            }
        }
    }


    private Map<Integer, DLinkedNode> cache = new HashMap<Integer, DLinkedNode>();
    private int size;
    private int capacity;
    private DLinkedNode head, tail;
    /* My Defined classes double linked list class */
    class DLinkedNode {
        int key;
        int value;
        DLinkedNode prev;
        DLinkedNode next;
    }

    private void addNode(DLinkedNode node) {
        //Since the new node is the most recently used node,
        //we add it after the pseudo head.
        //Pseudo Head and Tail of the DLL will remain the same.
        node.prev = head;
        node.next = head.next;

        //The former node after the pseudo head, now points to
        //the new node as its previous node.
        head.next.prev = node;
        head.next = node;
    }

    //Removing a node from the DLL
    // Do the same but in reverse
    private void removeNode(DLinkedNode node) {
        //The node before and after the deleted node
        DLinkedNode prev = node.prev;
        DLinkedNode next = node.next;

        //The previous node points to the node after
        //the removed node.
        prev.next = next;
        //The next node points to the node before
        //the removed node.
        next.prev = prev;
    }

    private void moveToHead(DLinkedNode node) {
        removeNode(node);
        addNode(node);
    }

    //Pops out the node before the pseudo tail
    //and returns the removed node.
    private DLinkedNode popTail() {
        DLinkedNode lastItem = tail.prev;
        lastItem.prev.next = tail;
        tail.prev = lastItem.prev;
        return lastItem;
    }

//    public LRUCache(int capacity) {
//        this.capacity = capacity;
//        size = 0;
//
//        tail = new DLinkedNode();
//        tail.key = -1;
//        tail.value = -1;
//        tail.next = null;
//
//        head = new DLinkedNode();
//        head.key = -1;
//        head.value = -1;
//        head.prev = null;
//        head.next = tail;
//
//        tail.prev = head;
//    }

    public int get(int key) {
        // if the cache contains the key, we return the item and move to head
        if (cache.containsKey(key)) {
            DLinkedNode item = cache.get(key);
            moveToHead(item);
            return item.value;
        } else {
            // If the cache does not contain the key
            return -1;
        }
    }

    public void put(int key, int value) {
        // if the cache contains the key, we update the key
        if (cache.containsKey(key)) {
            DLinkedNode item = cache.get(key);
            item.value = value;
            moveToHead(item);
        } else {
            // Add the new item
            DLinkedNode newItem = new DLinkedNode();
            newItem.key = key;
            newItem.value = value;
            cache.put(key, newItem);
            addNode(newItem);
            ++size;

            if (size > capacity) {
                // Get rid of the least recent item
                DLinkedNode lruItem = popTail();
                cache.remove(lruItem.key);
                --size;
            }
        }
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

    public class ListNode {
      int val;
      ListNode next;
      ListNode() {}
      ListNode(int val) { this.val = val; }
      ListNode(int val, ListNode next) { this.val = val; this.next = next; }
    }

    static class Trie {

        TrieObj rootNode;

        public class TrieObj {
            char val;
            HashMap<Character, TrieObj> dict = new HashMap<Character, TrieObj>();
            boolean isEndNode = false;
        }

        public Trie() {
            rootNode = new TrieObj();
        }

        public void insert(String word) {
            char[] wordCharArray = word.toCharArray();
            TrieObj currentNode = rootNode;

            for(int i = 0; i < wordCharArray.length; i++) {
                // If it does exist
                if (currentNode.dict.containsKey(wordCharArray[i])){

                    // Navigate to this node
                    currentNode = currentNode.dict.get(wordCharArray[i]);
                } else {
                    // if it doesn't exist we need to create it
                    TrieObj newNode = new TrieObj();
                    newNode.val = wordCharArray[i];

                    currentNode.dict.put(wordCharArray[i], newNode);
                    currentNode = newNode;
                }
            }

            // End the Current word
            currentNode.isEndNode = true;
        }

        public boolean search(String word) {
            char[] wordCharArray = word.toCharArray();
            TrieObj currentNode = rootNode;

            for (int i = 0; i < wordCharArray.length; i++) {
                if (!currentNode.dict.containsKey(wordCharArray[i])) {
                    return false;
                } else {
                    currentNode = currentNode.dict.get(wordCharArray[i]);
                }
            }

            if(currentNode.isEndNode) {
                return true;
            }

            return false;
        }

        public boolean startsWith(String prefix) {
            char[] wordCharArray = prefix.toCharArray();
            TrieObj currentNode = rootNode;

            for (int i = 0; i < wordCharArray.length; i++) {
                if (!currentNode.dict.containsKey(wordCharArray[i])) {
                    return false;
                }
            }

            return true;
        }
    }

}



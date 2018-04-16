我這次實做的是Problem 1.1~1.3 還有 Problem2 。

Problem1.1說明：檔名：prime.lsp
	檔案有寫main函式，首先執行程式後，使用者要輸入想要測試的數字，
然後就會印出他是不是質數。
實做方法是從 2 測試到 n/2，如果其中沒有數字可以被輸入的數字整除，就
表示他是質數。

輸入：17
輸出：17 is a prime number

Problem1.2說明：檔名：palindrome.lsp
	這個檔案同樣有main函式，執行程式後要輸入list，程式會判斷它
是不是list。 這裡我是用read讀進list，所以輸入格式要是 (m a d a m)
有左右括弧。
實做方法是用(equal l (reverse l))，看看本來的list是否跟reverse版本相同。

輸入：(cat dog bird bird dog cat)
輸出：(CAT DOG BIRD BIRD DOG CAT) is palindrome 
輸入範例二：()
輸出範例二：NIL is palindrome

Problem1.3說明：檔名：fib.lsp
	fib1是普通recursive版本，fib2是tail recursion版本。程式開始會先叫使用者
輸入想要的編號(編號從1開始)，代表數列的第幾個數字。例如輸入1 會印出0。
輸入2會印出1，輸入3會印出1，輸入4會印出2，輸入5會印出3......。
兩個實做方法不同之處在於fib2有用一個result來暫存答案，讓它不用一直呼叫很多層。
然後在main的部份，有先trace fib1 與trace fib2，所以後面呼叫的時候會印出
recursive的detail。此外最後也有印出兩個版本的答案供check。

輸入：5
輸出：
	(一連串recursive版本的資訊)
	Answer for orig: 3 
	-------------------------
	(一連串tail recursion版本的資訊)
	Answer for tail: 3


Problem2說明：檔名：Mergesort.lsp
	輸入格式皆與題目要求相同，先輸入有多少數字，再輸入各個數字。
	實做方法： defun mergesort是整個mergesort的主體，而merge的動作
在defun Mymerge裡面做處理，這裡是作merge的動作，比較兩個數列並排列。
輸入：
11
5 2 2 2 4 9 89 5 2 1 4 
輸出：
1 2 2 2 2 4 4 5 5 9 89
輸入範例2:
26
1 6 5 8 -5 6 9 89 54 54 3 7 4 10 52 65 98 54 52 46 5 62 6 100 99 0
輸出2：
-5 0 1 3 4 5 5 6 6 6 7 8 9 10 46 52 52 54 54 54 62 65 89 98 99 100
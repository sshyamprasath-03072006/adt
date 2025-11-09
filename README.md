# Complete C++ Guide for Data Structures & Algorithms

I'll cover all the major topics from your document with simple, easy-to-understand code examples.

---

## TOPIC 1: TIME COMPLEXITY & RECURSION

### 1. Default Parameters (calculateArea)

```cpp
#include <bits/stdc++.h>
using namespace std;

void calculateArea(int length, int breadth, char shape = 'r') {
    if (shape == 'r') {
        cout << "Area of rectangle: " << length * breadth << endl;
    } else if (shape == 't') {
        cout << "Area of triangle: " << length * breadth / 2 << endl;
    } else if (shape == 'c') {
        double area = 3.14 * length * length;
        cout << fixed << setprecision(2);
        cout << "Area of circle: " << area << endl;
    } else {
        cout << "Invalid shape!" << endl;
    }
}

int main() {
    int length, breadth;
    char shape;
    cin >> length >> breadth >> shape;
    calculateArea(length, breadth, shape);
    return 0;
}
```
**Walkthrough**: Default parameter `shape = 'r'` means if no shape is given, rectangle is assumed. Check shape and calculate accordingly.

---

### 2. Pass by Reference (Update Maximum)

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int a, b, newVal;
    cin >> a >> b >> newVal;
    
    cout << "Before: a = " << a << ", b = " << b << endl;
    
    int &maxRef = (a > b) ? a : b;  // Reference to max
    maxRef = newVal;
    
    cout << "After: a = " << a << ", b = " << b << endl;
    return 0;
}
```
**Walkthrough**: `&maxRef` creates a reference to the larger variable. Changing `maxRef` changes the original variable.

---

### 3. GCD using Recursion (Euclidean Algorithm)

```cpp
#include <bits/stdc++.h>
using namespace std;

int findGCD(int a, int b) {
    if (b == 0) return a;
    return findGCD(b, a % b);
}

int findLCM(int a, int b) {
    return (a * b) / findGCD(a, b);
}

int main() {
    int a, b;
    cin >> a >> b;
    cout << findGCD(a, b) << endl;
    cout << findLCM(a, b) << endl;
    return 0;
}
```
**Walkthrough**: GCD uses Euclidean algorithm recursively. LCM formula: `(a*b)/GCD(a,b)`.

---

### 4. Palindrome Check (Recursion)

```cpp
#include <bits/stdc++.h>
using namespace std;

int reverseNum(int n, int rev = 0) {
    if (n == 0) return rev;
    return reverseNum(n / 10, rev * 10 + n % 10);
}

int main() {
    int n;
    cin >> n;
    int rev = reverseNum(n);
    cout << "Reversed number: " << rev << endl;
    if (n == rev) {
        cout << "The number is a palindrome." << endl;
    } else {
        cout << "The number is not a palindrome." << endl;
    }
    return 0;
}
```
**Walkthrough**: Recursively build reversed number. Compare original with reversed.

---

### 5. Tower of Hanoi

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    int moves = (1 << n) - 1;  // 2^n - 1
    cout << moves << endl;
    return 0;
}
```
**Walkthrough**: Formula for minimum moves is `2^n - 1`. Use bit shift `(1 << n)` for `2^n`.

---

### 6. Sieve of Eratosthenes (NOT Sundaram)

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    
    vector<bool> isPrime(n + 1, true);
    isPrime[0] = isPrime[1] = false;
    
    for (int i = 2; i * i <= n; i++) {
        if (isPrime[i]) {
            for (int j = i * i; j <= n; j += i) {
                isPrime[j] = false;
            }
        }
    }
    
    for (int i = 2; i <= n; i++) {
        if (isPrime[i]) cout << i << " ";
    }
    return 0;
}
```
**Walkthrough**: Mark all multiples of primes as non-prime. Start from `i*i` for optimization.

---

## TOPIC 2: HEAPS

### 1. Min Heap (STL Priority Queue)

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    
    priority_queue<int, vector<int>, greater<int>> minHeap;
    int sum = 0, count = 0;
    
    for (int i = 0; i < n; i++) {
        int x;
        cin >> x;
        if (x > 0) {
            minHeap.push(x);
            sum += x;
            count++;
        }
    }
    
    if (count == 0) {
        cout << "No valid weight" << endl;
        return 0;
    }
    
    // Print heap
    priority_queue<int, vector<int>, greater<int>> temp = minHeap;
    while (!temp.empty()) {
        cout << temp.top() << " ";
        temp.pop();
    }
    cout << endl;
    
    cout << fixed << setprecision(2) << (double)sum / count << endl;
    return 0;
}
```
**Walkthrough**: Use `priority_queue` with `greater<int>` for min heap. Push positive values only.

---

### 2. Max Heap (Fibonacci Sequence)

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    
    priority_queue<int> maxHeap;
    int a = 2, b = 3;
    
    for (int i = 0; i < n; i++) {
        int curr = (i == 0) ? 2 : (i == 1) ? 3 : a + b;
        maxHeap.push(curr);
        
        cout << "Insert " << curr << ": ";
        priority_queue<int> temp = maxHeap;
        while (!temp.empty()) {
            cout << temp.top() << " ";
            temp.pop();
        }
        cout << endl;
        
        if (i >= 1) {
            a = b;
            b = curr;
        }
    }
    return 0;
}
```
**Walkthrough**: Default `priority_queue` is max heap. Generate Fibonacci and insert.

---

### 3. Heap Sort

```cpp
#include <bits/stdc++.h>
using namespace std;

void heapify(vector<int> &arr, int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    
    if (left < n && arr[left] > arr[largest]) largest = left;
    if (right < n && arr[right] > arr[largest]) largest = right;
    
    if (largest != i) {
        swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(vector<int> &arr) {
    int n = arr.size();
    for (int i = n / 2 - 1; i >= 0; i--) heapify(arr, n, i);
    for (int i = n - 1; i > 0; i--) {
        swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}

int main() {
    int n;
    cin >> n;
    vector<int> arr(n);
    for (int i = 0; i < n; i++) cin >> arr[i];
    
    heapSort(arr);
    
    for (int x : arr) cout << x << " ";
    return 0;
}
```
**Walkthrough**: Build max heap, then repeatedly extract max and heapify.

---

## TOPIC 3: GREEDY ALGORITHMS

### 1. Activity Selection

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    
    vector<int> start(n), finish(n);
    for (int i = 0; i < n; i++) cin >> start[i];
    for (int i = 0; i < n; i++) cin >> finish[i];
    
    vector<pair<int, int>> activities;
    for (int i = 0; i < n; i++) {
        activities.push_back({finish[i], i});
    }
    sort(activities.begin(), activities.end());
    
    cout << activities[0].second;
    int lastFinish = activities[0].first;
    
    for (int i = 1; i < n; i++) {
        int idx = activities[i].second;
        if (start[idx] >= lastFinish) {
            cout << " " << idx;
            lastFinish = activities[i].first;
        }
    }
    cout << endl;
    return 0;
}
```
**Walkthrough**: Sort by finish time. Greedily pick activities that don't overlap.

---

### 2. Fractional Knapsack

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int n;
    cin >> n;
    
    vector<int> weight(n), value(n);
    for (int i = 0; i < n; i++) cin >> weight[i];
    for (int i = 0; i < n; i++) cin >> value[i];
    
    int W;
    cin >> W;
    
    vector<pair<double, int>> ratio;
    for (int i = 0; i < n; i++) {
        ratio.push_back({(double)value[i] / weight[i], i});
    }
    sort(ratio.rbegin(), ratio.rend());
    
    double totalValue = 0;
    for (auto &p : ratio) {
        int idx = p.second;
        if (W >= weight[idx]) {
            totalValue += value[idx];
            cout << "Added object " << idx + 1 << " (Rs. " << value[idx] 
                 << ", " << weight[idx] << "Kg) completely in the bag. Space left: " 
                 << (W - weight[idx]) << "." << endl;
            W -= weight[idx];
        } else {
            int percent = (W * 100) / weight[idx];
            totalValue += (double)W * value[idx] / weight[idx];
            cout << "Added " << percent << "% (Rs." << value[idx] 
                 << ", " << weight[idx] << "Kg) of object " << idx + 1 << " in the bag." << endl;
            W = 0;
            break;
        }
    }
    
    cout << fixed << setprecision(2);
    cout << "Filled the bag with objects worth Rs. " << totalValue << "." << endl;
    return 0;
}
```
**Walkthrough**: Sort by value/weight ratio. Take items greedily, fractionally if needed.

---

### 3. Graph Coloring (Greedy)

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    int v;
    cin >> v;
    
    vector<vector<int>> graph(v, vector<int>(v));
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < v; j++) {
            cin >> graph[i][j];
        }
    }
    
    vector<int> color(v, 0);
    color[0] = 1;
    int maxColor = 1;
    
    for (int i = 1; i < v; i++) {
        set<int> usedColors;
        for (int j = 0; j < i; j++) {
            if (graph[i][j] == 1 && color[j] != 0) {
                usedColors.insert(color[j]);
            }
        }
        
        int c = 1;
        while (usedColors.count(c)) c++;
        color[i] = c;
        maxColor = max(maxColor, c);
    }
    
    cout << "Chromatic Number of the graph is: " << maxColor << endl;
    return 0;
}
```
**Walkthrough**: Assign smallest unused color to each vertex, avoiding neighbors' colors.

---

## TOPIC 4 & 5: STRING ALGORITHMS

### 1. Naive Pattern Matching

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    string text, pattern;
    getline(cin, text);
    getline(cin, pattern);
    
    int n = text.length(), m = pattern.length();
    
    for (int i = 0; i <= n - m; i++) {
        if (text.substr(i, m) == pattern) {
            cout << "Pattern found at index " << i << endl;
        }
    }
    return 0;
}
```
**Walkthrough**: Use `substr(i, m)` to extract substring of length m. Compare with pattern.

---

### 2. KMP Algorithm

```cpp
#include <bits/stdc++.h>
using namespace std;

vector<int> computeLPS(string pattern) {
    int m = pattern.length();
    vector<int> lps(m, 0);
    int len = 0, i = 1;
    
    while (i < m) {
        if (pattern[i] == pattern[len]) {
            lps[i++] = ++len;
        } else {
            if (len != 0) len = lps[len - 1];
            else lps[i++] = 0;
        }
    }
    return lps;
}

int main() {
    string text, pattern;
    getline(cin, text);
    getline(cin, pattern);
    
    vector<int> lps = computeLPS(pattern);
    int n = text.length(), m = pattern.length();
    int i = 0, j = 0;
    
    while (i < n) {
        if (text[i] == pattern[j]) {
            i++; j++;
        }
        
        if (j == m) {
            cout << "Found pattern at index " << i - j << endl;
            j = lps[j - 1];
        } else if (i < n && text[i] != pattern[j]) {
            if (j != 0) j = lps[j - 1];
            else i++;
        }
    }
    return 0;
}
```
**Walkthrough**: LPS array stores longest proper prefix which is also suffix. Use it to avoid re-checking.

---

### 3. Rabin-Karp Algorithm

```cpp
#include <bits/stdc++.h>
using namespace std;

int main() {
    string text, pattern;
    getline(cin, text);
    getline(cin, pattern);
    
    int n = text.length(), m = pattern.length();
    int d = 256, q = 101;
    int h = 1;
    
    for (int i = 0; i < m - 1; i++) h = (h * d) % q;
    
    int p = 0, t = 0;
    for (int i = 0; i < m; i++) {
        p = (d * p + pattern[i]) % q;
        t = (d * t + text[i]) % q;
    }
    
    int count = 0;
    for (int i = 0; i <= n - m; i++) {
        if (p == t) {
            bool match = true;
            for (int j = 0; j < m; j++) {
                if (text[i + j] != pattern[j]) {
                    match = false;
                    break;
                }
            }
            if (match) count++;
        }
        
        if (i < n - m) {
            t = (d * (t - text[i] * h) + text[i + m]) % q;
            if (t < 0) t += q;
        }
    }
    
    cout << count << endl;
    return 0;
}
```
**Walkthrough**: Rolling hash technique. Compare hash values first, then verify actual match.

---

### 4. Manacher's Algorithm (Longest Palindrome)

```cpp
#include <bits/stdc++.h>
using namespace std;

string longestPalindrome(string s) {
    string t = "#";
    for (char c : s) {
        t += c;
        t += '#';
    }
    
    int n = t.length();
    vector<int> p(n, 0);
    int center = 0, right = 0;
    
    for (int i = 0; i < n; i++) {
        int mirror = 2 * center - i;
        if (i < right) p[i] = min(right - i, p[mirror]);
        
        while (i + p[i] + 1 < n && i - p[i] - 1 >= 0 && 
               t[i + p[i] + 1] == t[i - p[i] - 1]) {
            p[i]++;
        }
        
        if (i + p[i] > right) {
            center = i;
            right = i + p[i];
        }
    }
    
    int maxLen = 0, centerIdx = 0;
    for (int i = 0; i < n; i++) {
        if (p[i] > maxLen) {
            maxLen = p[i];
            centerIdx = i;
        }
    }
    
    int start = (centerIdx - maxLen) / 2;
    return s.substr(start, maxLen);
}

int main() {
    string s;
    cin >> s;
    cout << longestPalindrome(s) << endl;
    return 0;
}
```
**Walkthrough**: Transform string with '#' separators. Use mirror property to expand around centers efficiently.

---

## TOPIC 6, 7, 8: BACKTRACKING

### 1. Rat in Maze

```cpp
#include <bits/stdc++.h>
using namespace std;

bool solve(vector<vector<int>> &maze, int x, int y, int n, 
           vector<vector<int>> &sol, string path, vector<string> &paths) {
    if (x == n - 1 && y == n - 1) {
        sol[x][y] = 1;
        paths.push_back(path);
        sol[x][y] = 0;
        return true;
    }
    
    if (x >= 0 && x < n && y >= 0 && y < n && maze[x][y] == 1 && sol[x][y] == 0) {
        sol[x][y] = 1;
        
        solve(maze, x + 1, y, n, sol, path + "D", paths);
        solve(maze, x, y - 1, n, sol, path + "L", paths);
        solve(maze, x, y + 1, n, sol, path + "R", paths);
        solve(maze, x - 1, y, n, sol, path + "U", paths);
        
        sol[x][y] = 0;
    }
    return false;
}

int main() {
    int n;
    cin >> n;
    
    vector<vector<int>> maze(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> maze[i][j];
        }
    }
    
    if (maze[0][0] == 0 || maze[n-1][n-1] == 0) {
        cout << -1 << endl;
        return 0;
    }
    
    vector<vector<int>> sol(n, vector<int>(n, 0));
    vector<string> paths;
    
    solve(maze, 0, 0, n, sol, "", paths);
    
    if (paths.empty()) {
        cout << -1 << endl;
    } else {
        sort(paths.begin(), paths.end());
        for (string p : paths) cout << p << " ";
    }
    return 0;
}
```
**Walkthrough**: Try all 4 directions (D, L, R, U). Backtrack if dead end. Store all valid paths.

---

### 2. N-Queens

```cpp
#include <bits/stdc++.h>
using namespace std;

bool isSafe(vector<vector<int>> &board, int row, int col, int n) {
    for (int i = 0; i < col; i++) {
        if (board[row][i]) return false;
    }
    for (int i = row, j = col; i >= 0 && j >= 0; i--, j--) {
        if (board[i][j]) return false;
    }
    for (int i = row, j = col; i < n && j >= 0; i++, j--) {
        if (board[i][j]) return false;
    }
    return true;
}

bool solve(vector<vector<int>> &board, int col, int n) {
    if (col >= n) return true;
    
    for (int i = 0; i < n; i++) {
        if (board[i][col] == 1) {
            return solve(board, col + 1, n);
        }
    }
    
    for (int i = 0; i < n; i++) {
        if (isSafe(board, i, col, n)) {
            board[i][col] = 1;
            if (solve(board, col + 1, n)) return true;
            board[i][col] = 0;
        }
    }
    return false;
}

int main() {
    int n, k;
    cin >> n >> k;
    
    vector<vector<int>> board(n, vector<int>(n, 0));
    for (int i = 0; i < k; i++) {
        int r, c;
        cin >> r >> c;
        board[r][c] = 1;
    }
    
    if (solve(board, 0, n)) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << board[i][j] << " ";
            }
            cout << endl;
        }
    } else {
        cout << "No solution found." << endl;
    }
    return 0;
}
```
**Walkthrough**: Place queens column by column. Check if safe (no row/diagonal conflicts). Backtrack if stuck.

---

### 3. Subset Sum (Partition)

```cpp
#include <bits/stdc++.h>
using namespace std;

void findSubsets(vector<int> &arr, int index, int target, vector<int> &curr, int &count) {
    if (target == 0) {
        count++;
        return;
    }
    if (index == arr.size() || target < 0) return;
    
    curr.push_back(arr[index]);
    findSubsets(arr, index + 1, target - arr[index], curr, count);
    curr.pop_back();
    
    findSubsets(arr, index + 1, target, curr, count);
}

int main() {
    int n;
    cin >> n;
    
    vector<int> arr(n);
    for (int i = 0; i < n; i++) cin >> arr[i];
    
    int target;
    cin >> target;
    
    vector<int> curr;
    int count = 0;
    findSubsets(arr, 0, target, curr, count);
    
    cout << count << endl;
    
    if (count % 2 == 0) {
        cout << count * 3 << endl;
    } else {
        cout << count + 7 << endl;
    }
    return 0;
}
```
**Walkthrough**: Try including/excluding each element. Count subsets with sum = target.

---

### 4. Sudoku Solver

```cpp
#include <bits/stdc++.h>
using namespace std;

bool isSafe(vector<vector<int>> &board, int row, int col, int num) {
    for (int x = 0; x < 9; x++) {
        if (board[row][x] == num || board[x][col] == num) return false;
    }
    
    int startRow = row - row % 3, startCol = col - col % 3;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (board[i + startRow][j + startCol] == num) return false;
        }
    }
    return true;
}

bool solveSudoku(vector<vector<int>> &board) {
    for (int row = 0; row < 9; row++) {
        for (int col = 0; col < 9; col++) {
            if (board[row][col] == 0) {
                for (int num = 1; num <= 9; num++) {
                    if (isSafe(board, row, col, num)) {
                        board[row][col] = num;
                        if (solveSudoku(board)) return true;
                        board[row][col] = 0;
                    }
                }
                return false;
            }
        }
    }
    return true;
}

int main() {
    vector<vector<int>> board(9, vector<int>(9));
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            cin >> board[i][j];
        }
    }
    
    if (solveSudoku(board)) {
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 9; j++) {
                cout << board[i][j] << " ";
            }
            cout << endl;
        }
    }
    return 0;
}
```
**Walkthrough**: Find empty cell. Try numbers 1-9. Check if safe (row, col, 3x3 box). Backtrack if no solution.

---

### 5. Graph Coloring

```cpp
#include <bits/stdc++.h>
using namespace std;

bool isSafe(int v, vector<vector<int>> &graph, vector<int> &color, int c) {
    for (int i = 0; i < graph.size(); i++) {
        if (graph[v][i] && c == color[i]) return false;
    }
    return true;
}

bool solve(vector<vector<int>> &graph, int m, vector<int> &color, int v) {
    if (v == graph.size()) return true;
    
    for (int c = 1; c <= m; c++) {
        if (isSafe(v, graph, color, c)) {
            color[v] = c;
            if (solve(graph, m, color, v + 1)) return true;
            color[v] = 0;
        }
    }
    return false;
}

int main() {
    int n;
    cin >> n;
    
    vector<vector<int>> graph(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> graph[i][j];
        }
    }
    
    int m;
    cin >> m;
    
    vector<int> color(n, 0);
    
    if (solve(graph, m, color, 0)) {
        cout << "Solution Exists:" << endl;
        for (int c : color) cout << c << " ";
    } else {
        cout << "Solution does not exist" << endl;
    }
    return 0;
}
```
**Walkthrough**: Try each color for vertex. Check if neighbors have same color. Backtrack if conflict.

---

## Key Tips for Exams:

1. **Use `#include <bits/stdc++.h>`** - includes everything
2. **STL is your friend**: `vector`, `priority_queue`, `set`, `map`
3. **String operations**: Use `substr()`, `find()`, `getline()`
4. **Recursion**: Always have base case first
5. **Backtracking**: Try option → Recurse → Undo (backtrack)
6. **Greedy**: Sort first, then make locally optimal choices
7. **Heaps**: Use `priority_queue` (default max, use `greater<int>` for min)

**Good luck with your exams!** 🎯

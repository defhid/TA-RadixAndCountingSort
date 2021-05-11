#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <chrono>

#define IS_DEMO false


using namespace std;

namespace io {

    template <class T>
    static void print_bites_of(T value, bool valuableOnly = true) {
        /** Output bites of value where left is the most significant bit */
        int bCount = sizeof(T) * 8;
        if (T(0) > T(-1)) {  // signed type
            if (value < T(0)) {
                cout << '-';
                value *= -1;
            }
            --bCount;
        }
        if (valuableOnly) {
            bool isValuable = false;
            for (int i = bCount - 1; i >= 0; --i) {
                if (value & (1 << i)) {
                    if (!isValuable)
                        isValuable = true;
                    cout << '1';
                }
                else {
                    if (isValuable) cout << '0';
                }
            }
        }
        else {
            for (int i = bCount - 1; i >= 0; --i) {
                cout << (value & (1 << i) ? '1' : '0');
            }
        }
    }

    template<class T>
    void print_vector(vector<T>& vec) {
        cout << "[";
        int end = vec.size() - 1;
        for (int i = 0; i <= end; ++i) {
            cout << +vec[i];  // '+' char fix
            if (i != end) cout << ", ";
        }
        cout << "]";
    }


    template<class T>
    void print_vector(const vector<T>& vec) {
        cout << "[";
        int end = vec.size() - 1;
        for (int i = 0; i <= end; ++i) {
            cout << +vec[i];
            if (i != end) cout << ", ";
        }
        cout << "]";
    }


    template<class T>
    void print_vector_binary(vector<T>& vec) {
        cout << "[";
        int end = vec.size() - 1;
        for (int i = 0; i <= end; ++i) {
            print_bites_of(vec[i], true);
            if (i != end) cout << ", ";
        }
        cout << "]";
    }


    inline void ask(const string& q) {
        if (cin.peek() == '\n') {
            cout << q;
        }
    }


    inline void clear_cin() {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }


    inline void make_fail(const string& m) {
        clear_cin();
        cout << endl << endl << m << endl;
    }
}



template<typename T>
class RadixSort {  // LSD
    /** T: char / short / int / long + unsigned */
public:
    static void numbers_binary(vector<T>& vec) {
#if IS_DEMO
        io::print_vector(vec);
        cout << endl;
        io::print_vector_binary(vec);
        cout << endl << endl;
#endif
        int vecSize = vec.size();
        if (vecSize < 2) return;

        int i, j;

        bool isSignedType = T(0) > T(-1);
        int bCount = sizeof(T) * 8;
        if (isSignedType) --bCount;  // minus sign bit (last)

        T* bitMasks = new T[bCount]();  // (1, 2, 4, 8, ...) or (1, 10, 100, 1000, ...)
        for (i = 0; i < bCount; ++i) {
            bitMasks[i] = pow(T(2), i);
        }

        T maxValue = *max_element(vec.begin(), vec.end());
        while (bCount > 1 && !(maxValue & bitMasks[bCount - 1])) {  // find max required mask (bCount)
            --bCount;
        }

        int lastNegative = _separate(vec, vecSize);
        int negCount = lastNegative + 1;
        int posCount = vecSize - negCount;

        int zerosCount, onesCount;
        auto tmp = vector<T>((posCount > negCount ? posCount : negCount), 0);

        for (i = 0; i < bCount; ++i) {
            T mask = bitMasks[i];
            zerosCount = 0;  // use as the current filling index
            onesCount = 0;
            for (j = 0; j < negCount; ++j) {
                if (vec[j] & mask) {  // 1
                    tmp[onesCount++] = vec[j];
                }
                else {  // 0
                    vec[zerosCount++] = vec[j];
                }
            }
            for (j = 0; j < onesCount; ++j) {
                vec[zerosCount++] = tmp[j];
            }
            onesCount = 0;
            for (j = negCount; j < vecSize; ++j) {
                if (vec[j] & mask) {  // 1
                    tmp[onesCount++] = vec[j];
                }
                else {  // 0
                    vec[zerosCount++] = vec[j];
                }
            }
            for (j = 0; j < onesCount; ++j) {
                vec[zerosCount++] = tmp[j];
            }
#if IS_DEMO
            cout << "Bit: ";
            io::print_bites_of(mask, false);
            cout << endl;
            io::print_vector(vec);
            cout << endl;
            io::print_vector_binary(vec);
            cout << endl << endl;
#endif
        }
        delete[] bitMasks;
    }

    static void numbers_by_remainder(vector<T>& vec) {
        int vecSize = vec.size();
        if (vecSize < 2) return;

        int i, j, d, tmp;
        int count[10] = { 0 };
        auto digits = vector<T>(vecSize, 0);
        auto sorted = vector<T>(vecSize, 0);

        int lastNegative = _separate(vec, vecSize);
        int negCount = lastNegative + 1;

        int maxLength = to_string(*max_element(vec.begin(), vec.end())).length();

        for (i = 0; i < maxLength; ++i) {
            // NEGATIVE
            for (j = 0; j < negCount; ++j) {
                d = -int((vec[j] / T(pow(10, i))) % 10);
                ++count[d];
                digits[j] = d;
            }
            tmp = count[9];
            count[9] = 0;
            for (j = 8; j >= 0; --j) {  // make indexes from counters (offset)
                swap(tmp, count[j]);
                count[j] += count[j + 1];
            }
            for (j = 0; j < negCount; ++j) {
                sorted[count[digits[j]]++] = vec[j];
            }
            for (j = 0; j < 10; ++j) {
                count[j] = negCount;
            }
            // POSITIVE
            for (j = negCount; j < vecSize; ++j) {
                d = (vec[j] / T(pow(10, i))) % 10;
                ++count[d];
                digits[j] = d;
            }
            tmp = count[0];
            count[0] = negCount;
            for (j = 1; j < 10; ++j) {
                swap(tmp, count[j]);
                count[j] += count[j - 1] - negCount;
            }
            for (j = negCount; j < vecSize; ++j) {
                sorted[count[digits[j]]++] = vec[j];
            }
            for (j = 0; j < 10; ++j) {
                count[j] = 0;
            }
            swap(vec, sorted);
        }
    }

    static void numbers_as_strings(vector<T>& vec) {
        int vecSize = vec.size();
        if (vecSize < 2) return;

        int i, j, tmp;
        int count[10] = { 0 };
        auto vecStrings = vector<string>();
        vecStrings.reserve(vecSize);

        int lastNegative = _separate(vec, vecSize);
        int negCount = lastNegative + 1;

        // fill strings vector and find max length among all strings
        size_t maxLength = 0;
        size_t currLength;
        string currStr;
        for (i = 0; i < vecSize; ++i) {
            currStr = to_string((i < negCount ? vec[i] * T(-1) : vec[i]));
            vecStrings.push_back(currStr);
            currLength = currStr.length();
            if (currLength > maxLength) {
                maxLength = currLength;
            }
        }

        // add '0' to make the same length
        for (i = 0; i < vecSize; ++i) {
            currLength = vecStrings[i].length();
            if (currLength < maxLength) {
                vecStrings[i].insert(0, maxLength - currLength, '0');
            }
        }

        auto sortedStrings = vector<string>(vecStrings);

        for (i = int(maxLength) - 1; i >= 0; --i) {
            // NEGATIVE
            for (j = 0; j < negCount; ++j) {
                ++count[_get_digit(vecStrings[j][i])];
            }
            tmp = count[9];
            count[9] = 0;
            for (j = 8; j >= 0; --j) {  // make indexes from counters (offset)
                swap(tmp, count[j]);
                count[j] += count[j + 1];
            }
            for (j = 0; j < negCount; ++j) {
                sortedStrings[count[_get_digit(vecStrings[j][i])]++] = vecStrings[j];
            }
            for (j = 0; j < 10; ++j) {
                count[j] = negCount;
            }
            // POSITIVE
            for (j = negCount; j < vecSize; ++j) {
                ++count[_get_digit(vecStrings[j][i])];
            }
            tmp = count[0];
            count[0] = negCount;
            for (j = 1; j < 10; ++j) {
                swap(tmp, count[j]);
                count[j] += count[j - 1] - negCount;
            }
            for (j = negCount; j < vecSize; ++j) {
                sortedStrings[count[_get_digit(vecStrings[j][i])]++] = vecStrings[j];
            }
            for (j = 0; j < 10; ++j) {
                count[j] = 0;
            }
            swap(vecStrings, sortedStrings);
        }

        if (std::is_same<T, long long>::value) {
            for (j = 0; j < negCount; ++j) {
                vec[j] = -stoll(vecStrings[j]);
            }
            for (j = negCount; j < vecSize; ++j) {
                vec[j] = stoll(vecStrings[j]);
            }
        }
        else if (std::is_same<T, unsigned long long>::value) {
            for (j = negCount; j < vecSize; ++j) {
                vec[j] = stoull(vecStrings[j]);
            }
        }
        else {
            for (j = 0; j < negCount; ++j) {
                vec[j] = -stoi(vecStrings[j]);
            }
            for (j = negCount; j < vecSize; ++j) {
                vec[j] = stoi(vecStrings[j]);
            }
        }
    }

private:
    static int _separate(vector<T>& vec, int vecSize) {
        /** Moves all negative numbers to the left side of vector.
         *  Returns: index of the last negative number (or -1). */
        int n = vecSize - 1;
        int i, j;
        bool found;
        for (i = 0; i < n; ++i) {
            if (vec[i] > 0) {
                found = false;
                for (j = i + 1; j < vecSize; ++j) {
                    if (vec[j] < 0) {
                        swap(vec[i], vec[j]);
                        found = true;
                    }
                }
                if (!found) return i - 1;
            }
        }
        cout << "!!! ";
        io::print_vector(vec);
        cout << " !!!" << endl;
        return i;
    }

    inline static int _get_digit(char c) {
        /** Returns: integer from 0 to 9 */
        return int(c - 48);  // '0' = 48, '1' = 49, ...
    }
};



class CountingSort {
public:
    template <typename T>
    static void v1(vector<T>& vec) {
        /** T: char / short / int / long + unsigned */
        int vecSize = vec.size();
        if (vecSize < 2) return;

        T maxValue = vec[0];
        T minValue = vec[0];
        for (int i = 0; i < vecSize; ++i) {
            if (vec[i] > maxValue) {
                maxValue = vec[i];
            }
            if (vec[i] < minValue) {
                minValue = vec[i];
            }
        }

        unsigned long long tmpRange = maxValue - minValue + 1;
        if (tmpRange > 70'000ull) {
            throw bad_exception();
        }
        int range = int(tmpRange);
        auto count = new T[range]{ T(0) };

        for (T value : vec) {
            count[value - minValue] += 1;
        }
        int last = 0;

#if IS_DEMO
        io::print_vector(vec);
        cout << endl << endl;
        vec.clear();
        vec.reserve(vecSize);
        for (int i = 0; i < range; ++i) {
            if (count[i] > 0) {
                for (; count[i] > 0; --count[i]) {
                    vec.push_back(minValue + i);
                }
                io::print_vector(vec);
                cout << endl << endl;
            }
        }
#else
        for (int i = 0; i < range; ++i) {
            if (count[i] > 0) {
                for (; count[i] > 0; --count[i]) {
                    vec[last++] = minValue + i;
                }
            }
        }
#endif
    }

    template <class T>
    static void v2(vector<T>& vec) {
        /** T: any type supported '>', '<' operators */
        int vecSize = vec.size();
        if (vecSize < 2) return;

        int i, j, pos;
        auto result = vector<T>(vecSize, 0);

#if IS_DEMO
        io::print_vector(vec);
        cout << endl << endl;
#endif
        for (i = 0; i < vecSize; ++i) {
            pos = 0;
            for (j = 0; j < i; ++j) {  // before vec[i] (with same values counting)
                if (vec[j] <= vec[i]) ++pos;
            }
            for (j = i + 1; j < vec.size(); ++j) {  // after vec[i]
                if (vec[j] < vec[i]) ++pos;
            }
            result[pos] = vec[i];
#if IS_DEMO
            io::print_vector(result);
            cout << endl << endl;
#endif
        }
        vec = result;
    }
};



template <class T>
class QuickSort {
    /** T: any type supported '>', '<' operators */
public:
    static void recursive(vector<T>& vec) {
        _recursive(vec, 0, vec.size() - 1);
    }

    static void cpp(vector<T>& vec) {
        sort(vec.begin(), vec.end());
    }

private:
    static void _recursive(vector<T>& vec, int left, int right) {
        if (left < right) {
            int pivotIndex = left + (right - left) / 2;
            T pivotValue = vec[pivotIndex];
            int i = left, j = right;
            while (i <= j) {
                while (vec[i] < pivotValue) {
                    ++i;
                }
                while (vec[j] > pivotValue) {
                    --j;
                }
                if (i <= j) {
                    swap(vec[i], vec[j]);
                    ++i;
                    --j;
                }
            }
            _recursive(vec, left, i - 1);
            _recursive(vec, i, right);
        }
    }
};


template <typename T>
class SortComparer {
public:
    struct Result {
    public:
        string name;
        vector<T> vec;
        unsigned long long time_ns;
        Result(string& sortName, vector<T> vecSorted, SortComparer* comp) : name(sortName), vec(std::move(vecSorted)) {
            time_ns = (comp->_endTime - comp->_startTime).count() * 1'000'000'000;
            comp->_results.push_back(*this);
        }
    };

private:
    vector<Result> _results;
    vector<T> _vec;  // origin
    chrono::time_point<chrono::steady_clock, chrono::duration<double>> _startTime;
    chrono::time_point<chrono::steady_clock, chrono::duration<double>> _endTime;

public:
    explicit SortComparer(vector<T> vecOrigin) : _vec(std::move(vecOrigin)) { }

    Result& execute(string sortName, void (sortFunc)(vector<T>&)) {  // make copy and sort it
        auto vec = vector<T>(_vec);
        _startTime = chrono::steady_clock::now();
        sortFunc(vec);
        _endTime = chrono::steady_clock::now();
        Result(sortName, std::move(vec), this);
        return _results.back();
    }

    vector<T>& getOriginVec() { return _vec; }
    vector<Result>& getResults() { return _results; }

    void print_results(bool printVec = false) {
        if (printVec) {
            cout << "SOURCE ARRAY:  ";
            io::print_vector(_vec);
            cout << endl << endl;
        }

        int maxSortNameLength = 0;
        int maxTime = 0;
        for (const Result& result : _results) {
            if (result.name.length() > maxSortNameLength) {
                maxSortNameLength = result.name.length();
            }
            if (result.time_ns > maxTime) {
                maxTime = result.time_ns;
            }
        }

        int maxTimeStringLength = to_string(maxTime).length();
        if (maxTimeStringLength % 3 == 0)
            maxTimeStringLength += int(maxTimeStringLength / 3) - 1;
        else
            maxTimeStringLength += int(maxTimeStringLength / 3);

        string timeString;
        int i, spaces;
        for (const Result& result : _results) {
            timeString = to_string(result.time_ns);
            cout << result.name << ":  ";
            spaces = (maxSortNameLength - result.name.length()) +
                     (maxTimeStringLength - timeString.length()) - int(timeString.length() / 3);
            cout << string(spaces, ' ');
            for (i = 0; i < timeString.length(); ++i) {
                if ((timeString.length() - i) % 3 == 0) {
                    cout << ' ';
                }
                cout << timeString[i];
            }
            cout << " ns";
            if (printVec) {
                cout << "  ";
                io::print_vector(result.vec);
                cout << endl << endl;
            }
            else cout << endl;
        }
    }
};


template <typename T>
vector<T> launch_compare(vector<T>& vec, bool fullResult) {
    auto comp = SortComparer<T>(std::move(vec));

    comp.execute("Radix binary", RadixSort<T>::numbers_binary);

    comp.execute("Radix by remainder", RadixSort<T>::numbers_by_remainder);

    comp.execute("Radix by strings", RadixSort<T>::numbers_as_strings);

    try { comp.execute("Counting v1.0", CountingSort::v1<T>); }
    catch (bad_exception&) {}

    comp.execute("Counting v2.0", CountingSort::v2<T>);

    comp.execute("Quick", QuickSort<T>::recursive);

    auto intrasortResult = comp.execute("Introsort (STL algorithm)", QuickSort<T>::cpp);

    comp.print_results(fullResult);
    return std::move(intrasortResult.vec);
}


template <typename T>
void launch_demo(vector<T>& vec1) {
    auto vec2 = vector<T>(vec1);
    auto vec3 = vector<T>(vec1);

    cout << "\n\n---RADIX BINARY SORT---\n\n";
    RadixSort<T>::numbers_binary(vec1);

    cout << "\n\n---COUNTING SORT v1.0---\n\n";
    try { CountingSort::v1(vec2); }
    catch (bad_exception&) { cout << "[Not suitable]" << endl; }

    cout << "\n\n---COUNTING SORT v2.0---\n\n";
    CountingSort::v2(vec3);

    cout << endl;
}


template <typename T>
void launch() {
#if !IS_DEMO
    io::ask("Select comparison mode (short/full): ");
    string mode;
    cin >> mode;

    bool isShort = (mode == "short" || mode == "s");
    bool isFull = (mode == "full" || mode == "f");

    if (!(isShort || isFull)) return io::make_fail("Incorrect mode!");
#endif

    io::ask("Input array size (1 <= size <= 100000): ");
    int n;
    cin >> n;

    if (cin.fail()) return io::make_fail("Incorrect array size!");
    if (n < 1 || n > 100'000) return io::make_fail("Incorrect array size! (" + to_string(n) + ')');

    io::ask("Input array nums:\n");
    auto vec = vector<T>();
    vec.reserve(n);

    if (std::is_same<T, char>::value || std::is_same<T, unsigned char>::value) {
        int num;
        for (int i = 0; i < n; ++i) {
            cin >> num;
            if (cin.fail()) {
                string val;
                cin >> val;
                return io::make_fail("Incorrect num! ('" + val + "')");
            }
            vec.push_back(T(num));
        }
    }
    else {
        T num;
        for (int i = 0; i < n; ++i) {
            cin >> num;
            if (cin.fail()) {
                string val;
                cin >> val;
                return io::make_fail("Incorrect num! ('" + val + "')");
            }
            vec.push_back(num);
        }
    }

#if IS_DEMO
    launch_demo(vec);
#else
    cout << endl << endl;
    if (!isFull) cout << "SOURCE ARRAY" << endl << endl;
    auto vecSorted = launch_compare<T>(vec, isFull);

    cout << endl;
    if (!isFull) cout << "SORTED ARRAY" << endl << endl;
    launch_compare<T>(vecSorted, isFull);
#endif
}


int main() {
#if IS_DEMO
    cout << " ----------------------------------------\n";
    cout << "|  RADIX & BINARY SORTING DEMONSTRATION  |\n";
    cout << " ---------------------------------------- \n\n";
#else
    cout << " ------------------------------\n";
    cout << "|  SORTING METHODS COMPARISON  |\n";
    cout << " ------------------------------ \n\n";
#endif
    cout << "Select nums type (char/short/int/long/uchar/ushort/uint/ulong): ";
    string type;
    cin >> type;
    if (type == "char" || type == "c") {
        launch<char>();
    }
    else if (type == "short" || type == "s") {
        launch<short>();
    }
    else if (type == "int" || type == "i") {
        launch<int>();
    }
    else if (type == "long" || type == "l") {
        launch<long long>();
    }
    else if (type == "uchar" || type == "uc") {
        launch<unsigned char>();
    }
    else if (type == "ushort" || type == "us") {
        launch<unsigned short>();
    }
    else if (type == "uint" || type == "ui") {
        launch<unsigned int>();
    }
    else if (type == "ulong" || type == "ul") {
        launch<unsigned long long>();
    }
    else {
        io::make_fail("Incorrect nums type! ('" + type + "')");
    }
    cout << endl << "Input anything to exit...";
    io::clear_cin();
    cin.get();
    return 0;
}
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <iterator>
#include "compact_vector.hpp"
#include <math.h>
#include <ctime>
#include <bit>
#include <bitset>
#include <cstdint>
#include <typeinfo>
#include <chrono>

class Bitvector {
public:
    Bitvector();
    Bitvector(std::string bitvector_string);
    bool save(std::string fname);
    bool load(std::string fname);
    int get_length();
    compact::vector<uint64_t, 1> get_bitvector();
    std::string get_b_string();
    void set_bitvector(std::string bitvector_string);
private:
    int length;
    compact::vector<uint64_t, 1> bitvector;
    std::string b_string;
};

Bitvector::Bitvector() {
    compact::vector<uint64_t, 1> u(0);
    std::swap(bitvector, u);
    b_string = "";
    length = 0;
}

Bitvector::Bitvector(std::string bitvector_string) {
    int bitvector_string_length = bitvector_string.length();
    // Here want to append enough 0's at end so that length is an even power of 2.
    // Note that this will increase the length of the bitvector by at most 4 (that is, by a constant
    // multiple)
    int new_length;
    new_length = ceil(log2(bitvector_string_length)) + (int (ceil(log2(bitvector_string_length))) % 2);
    new_length = pow(2, new_length);
    int num_zeros_to_add = new_length - bitvector_string_length;
    for (int i = 0; i < num_zeros_to_add; i++) {
        bitvector_string += '0';
    }
    compact::vector<uint64_t, 1> u(new_length);
    std::swap(bitvector, u);
    length = new_length;
    b_string = bitvector_string;
    for (int i = 0; i < new_length; i++) {
        bitvector[i] = bitvector_string[i];
    }
}

int Bitvector::get_length() {
    return length;
}

compact::vector<uint64_t, 1> Bitvector::get_bitvector() {
    return bitvector;
}

std::string Bitvector::get_b_string() {
    return b_string;
}

void Bitvector::set_bitvector(std::string bitvector_string) {
    int bitvector_string_length = bitvector_string.length();
    // Here want to append enough 0's at end so that length is an even power of 2.
    // Note that this will increase the length of the bitvector by at most 4 (that is, by a constant
    // multiple)
    int new_length;
    new_length = ceil(log2(bitvector_string_length)) + (int (ceil(log2(bitvector_string_length))) % 2);
    new_length = pow(2, new_length);
    int num_zeros_to_add = new_length - bitvector_string_length;
    for (int i = 0; i < num_zeros_to_add; i++) {
        bitvector_string += '0';
    }
    compact::vector<uint64_t, 1> u(new_length);
    std::swap(bitvector, u);
    length = new_length;
    b_string = bitvector_string;
    for (int i = 0; i < new_length; i++) {
        bitvector[i] = bitvector_string[i];
    }
}

bool Bitvector::save(std::string fname) {
    std::ofstream outputFile(fname);
    outputFile << b_string << std::endl;
    outputFile << length << std::endl;
    outputFile.close();
    return true;
}

bool Bitvector::load(std::string fname) {
    std::ifstream inputFile;
    inputFile.open(fname);
    if (!inputFile) {
        std::cout << "The file could not be opened.\n";
        return false;
    }

    getline(inputFile, b_string);
    std::string l;
    getline(inputFile, l);
    length = std::stoi(l);

    compact::vector<uint64_t, 1> u(length);
    std::swap(bitvector, u);
    for (int i = 0; i < length; i++) {
        bitvector[i] = b_string[i];
    }

    return true;
}

class Rank {
public:
    Rank();
    Rank(Bitvector * b);
    uint64_t rank1(uint64_t i);
    uint64_t overhead();
    bool save(std::string fname);
    bool load(std::string fname);
    int get_n();
private:
    Bitvector bitvector;
    std::vector<int> popcount_of_chunks;
    std::vector<std::vector<int>> popcount_of_subchunks;
    uint64_t ohead;
    // nc = num_chunks, ns = num_subchunks per chunk, lc = length of chunk, ls = length of subchunk
    // n is length of bitvector
    // ns_l is the number of subchunks in the last chunk
    // lc_l is the length of the last chunk
    int nc, ns, lc, ls, n, ns_l, lc_l;
};

// Default constructor
Rank::Rank() {
    popcount_of_chunks = {};
    popcount_of_subchunks = {{}};
    ohead = 0;
    nc = 0;
    ns = 0;
    lc = 0;
    ls = 0;
    n = 0;
    ns_l = 0;
    lc_l = 0;
}

// Constructor given bitvector
Rank::Rank(Bitvector * b) {
    bitvector = *b;
    // Here calculate chunks, subchunks, and ohead
    n = bitvector.get_length();
    lc = pow(log2(n),2);
    nc = ceil(static_cast<float>(n) / lc);
    ns = 2*log2(n);
    ls = log2(n)/2;
    int chunk_popcount = 0;
    int subchunk_popcount = 0;
    for (int i = 0; i < nc-1; i++) {
        popcount_of_chunks.push_back(chunk_popcount);
        subchunk_popcount = 0;
        std::vector<int> subchunk_popcount_vector;
        for (int j = 0; j < ns; j++) {
            subchunk_popcount_vector.push_back(subchunk_popcount);
            // popcount of subchunk
            subchunk_popcount += std::popcount(bitvector.get_bitvector().get_int(i*lc+j*ls,ls));
        }
        // add up popcounts of subchunks to get popcount of chunk and input that into correct part of chunk vector
        chunk_popcount += subchunk_popcount;
        popcount_of_subchunks.push_back(subchunk_popcount_vector);
    }
    popcount_of_chunks.push_back(chunk_popcount);
    subchunk_popcount = 0;
    std::vector<int> subchunk_popcount_vector;
    lc_l = n % lc;
    ns_l = lc_l/ls;
    for (int j = 0; j < ns_l; j++) {
        subchunk_popcount_vector.push_back(subchunk_popcount);
        subchunk_popcount += std::popcount(bitvector.get_bitvector().get_int((nc-1)*lc+j*ls,ls));
    }
    // add up popcounts of subchunks to get popcount of chunk and input that into correct part of chunk vector
    chunk_popcount += subchunk_popcount;
    popcount_of_subchunks.push_back(subchunk_popcount_vector);

    // calculate ohead
    // 32 is the number of bits taken up by an int in memory
    ohead = 32*popcount_of_chunks.size()*(1 + popcount_of_subchunks[0].size());
}

uint64_t Rank::overhead() {
    return ohead;
}

uint64_t Rank::rank1(uint64_t i) {
    // Here calculate rank
    uint64_t rank;

    // this is the case where an element is queried which is beyond the length of the bitvector
    assert(i < n);

    int q1 = floor(i/lc);
    int r = i % lc;
    int q2 = floor(r/ls);
    int cum_rank = 0;
    cum_rank = popcount_of_chunks[q1];
    int rel_cum_rank = 0;
    rel_cum_rank = popcount_of_subchunks[q1][q2];
    int rank_within_subchunk = 0;
    rank_within_subchunk = std::popcount(bitvector.get_bitvector().get_int(q1*lc+q2*ls,i-(q1*lc+q2*ls)));

    rank = cum_rank + rel_cum_rank + rank_within_subchunk;

    return rank;
}

bool Rank::save(std::string fname) {
    std::ofstream outputFile(fname);
    outputFile << n << std::endl << nc << std::endl << lc << std::endl << ns << std::endl << ls << std::endl;
    outputFile << ns_l << std::endl << lc_l << std::endl;
    outputFile << ohead << std::endl;
    outputFile << bitvector.get_b_string() << std::endl;
    // Now still need to save popcount of chunks, popcount of subchunks
    for (int i = 0; i < nc; i++) {
        outputFile << popcount_of_chunks[i] << std::endl;
    }
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < ns; j++) {
            outputFile << popcount_of_subchunks[i][j] << std::endl;
        }
    }
    outputFile.close();
    return true;
}

bool Rank::load(std::string fname) {
    std::ifstream inputFile;
    inputFile.open(fname);
    if (!inputFile) {
        std::cout << "The file could not be opened.\n";
        return false;
    }
    std::string n_str, nc_str, lc_str, ns_str, ls_str, ns_l_str, lc_l_str, ohead_str;
    getline(inputFile, n_str);
    n = std::stoi(n_str);
    getline(inputFile, nc_str);
    nc = std::stoi(nc_str);
    getline(inputFile, lc_str);
    lc = std::stoi(lc_str);
    getline(inputFile, ns_str);
    ns = std::stoi(ns_str);
    getline(inputFile, ls_str);
    ls = std::stoi(ls_str);
    getline(inputFile, ns_l_str);
    ns_l = std::stoi(ns_l_str);
    getline(inputFile, lc_l_str);
    lc_l = std::stoi(lc_l_str);
    getline(inputFile, ohead_str);
    ohead = std::stoi(ohead_str);

    std::string bitvector_str;
    getline(inputFile, bitvector_str);
    bitvector.set_bitvector(bitvector_str);

    // now popcount of chunks, popcount of subchunks
    int chunk_popcount;
    std::string chunk_popcount_str;
    for (int i = 0; i < nc; i++) {
        getline(inputFile, chunk_popcount_str);
        chunk_popcount = std::stoi(chunk_popcount_str);
        popcount_of_chunks.push_back(chunk_popcount);
    }

    int subchunk_popcount;
    std::string subchunk_popcount_str;
    for (int i = 0; i < nc; i++) {
        std::vector<int> subchunk_popcount_vector;
        for (int j = 0; j < ns; j++) {
            getline(inputFile, subchunk_popcount_str);
            subchunk_popcount = std::stoi(subchunk_popcount_str);
            subchunk_popcount_vector.push_back(subchunk_popcount);

        }
        popcount_of_subchunks.push_back(subchunk_popcount_vector);
    }

    return true;

}

int Rank::get_n() {
    return n;
}

class Select {
public:
    Select();
    Select(Rank * r);
    uint64_t select1(uint64_t i);
    uint64_t overhead();
    bool save(std::string fname);
    bool load(std::string fname);
private:
    uint64_t ohead;
    Rank rank_obj;
};

Select::Select() {
    rank_obj = Rank();
    ohead = 0;
}

Select::Select(Rank * r) {
    rank_obj = *r;
    ohead = rank_obj.overhead();
}

uint64_t Select::select1(uint64_t i) {
    uint64_t j;
    int n = rank_obj.get_n();
    assert(n > 0);
    // This checks the condition to see if there is no element of the bitvector
    // whose rank is the given argument
    if (rank_obj.rank1(n-1) < i) {
        std::cout << "No answer";
        return -1;
    }
    int l, c, r;
    l = 0;
    r = n-1;
    int found_answer = 0;
    while (found_answer == 0) {
        c = floor((l+r)/2);
        if (rank_obj.rank1(c) == i) {
            found_answer = 1;
        }
        else if (rank_obj.rank1(c) > i) {
            r = c;
        }
        else {
            l = c;
        }
    }

    j = c;
    while (j >= 0 && rank_obj.rank1(j) == i) {
        j--;
    }
    j++;

    return j;
}

uint64_t Select::overhead() {
    return ohead;
}

bool Select::save(std::string fname) {
    rank_obj.save(fname);
    return true;
}

bool Select::load(std::string fname) {
    rank_obj.load(fname);
    ohead = rank_obj.overhead();
    return true;
}

class Sparse {
public:
    Sparse();
    Sparse(uint64_t size);
    void append(std::string elem, uint64_t pos);
    void finalize();
    bool get_at_rank(uint64_t r, std::string& elem);
    bool get_at_index(uint64_t r, std::string& elem);
    uint64_t get_index_of(uint64_t r);
    uint64_t num_elem_at(uint64_t r);
    uint64_t size();
    uint64_t num_elem();
    bool save(std::string fname);
    bool load(std::string fname);
private:
    Bitvector bitvector;
    std::string bitvector_str;
    std::vector<std::string> values;
    Rank rank_obj;
    Select select_obj;
    bool finalized;
    uint64_t sz;
};

Sparse::Sparse(uint64_t size) {
    sz = size;
    bitvector_str = std::string(sz, '0');
    bitvector = Bitvector(bitvector_str);
}

void Sparse::append(std::string elem, uint64_t pos) {
    if (finalized == false) {
        values.push_back(elem);
        bitvector_str[pos] = '1';
        bitvector.set_bitvector(bitvector_str);
    }
}

void Sparse::finalize() {
    finalized = true;
    rank_obj = Rank(&bitvector);
    select_obj = Select(&rank_obj);
}

bool Sparse::get_at_rank(uint64_t r, std::string&  elem) {
    if (values.size() >= r) {
        elem = values[r];
        return true;
    }
    else {
        return false;
    }
}

bool Sparse::get_at_index(uint64_t r, std::string& elem) {
    if (bitvector.get_length() > r+1) {
        if (bitvector.get_bitvector()[r] == 1) {
            elem = values[rank_obj.rank1(r)];
            return true;
        } else {
            return false;
        }
    }
    else {
        return false;
    }
}

uint64_t Sparse::get_index_of(uint64_t r) {
    if (rank_obj.rank1(rank_obj.get_n()-1) >= r) {
        return select_obj.select1(r) - 1;
    }
    else if (bitvector.get_bitvector()[bitvector.get_length()-1] == 1 && rank_obj.rank1(rank_obj.get_n()-1) == r-1) {
        return bitvector.get_length()-1;
    }
    else {
        return -1;
    }
}

uint64_t Sparse::num_elem_at(uint64_t r) {
    if (bitvector.get_bitvector()[r] == 1) {
        return rank_obj.rank1(r) + 1;
    }
    else {
        return rank_obj.rank1(r);
    }
}

uint64_t Sparse::size() {
    return sz;
}

uint64_t Sparse::num_elem() {
    return values.size();
}

bool Sparse::save(std::string fname) {
    std::ofstream outputFile(fname);
    outputFile << bitvector_str << std::endl;
    for (int i = 0; i < values.size(); i++) {
        outputFile << values[i] << std::endl;
    }
    outputFile.close();
    return true;
}
bool Sparse::load(std::string fname) {
    std::ifstream inputFile;
    inputFile.open(fname);
    if (!inputFile) {
        std::cout << "The file could not be opened.\n";
        return false;
    }
    getline(inputFile, bitvector_str);
    bitvector = Bitvector(bitvector_str);
    finalized = true;
    rank_obj = Rank(&bitvector);
    select_obj = Select(&rank_obj);
    int bitvector_length = bitvector.get_length();
    sz = this->num_elem_at(bitvector_length-1);
    for (int i = 0; i < sz; i++) {
        std::string value;
        getline(inputFile, value);
        values.push_back(value);
    }
    sz = bitvector.get_length();
    return true;
}

int main() {

    std::string b_str_1000_1, b_str_1000_5, b_str_1000_10;
    std::string b_str_10000_1, b_str_10000_5, b_str_10000_10;
    std::string b_str_100000_1, b_str_100000_5, b_str_100000_10;
    std::string b_str_1000000_1, b_str_1000000_5, b_str_1000000_10;

    b_str_1000_1 = b_str_1000_5 = b_str_1000_10 = std::string(1000, '0');
    b_str_10000_1 = b_str_10000_5 = b_str_10000_10 = std::string(10000, '0');
    b_str_100000_1 = b_str_100000_5 = b_str_100000_10 = std::string(100000, '0');
    b_str_1000000_1 = b_str_1000000_5 = b_str_1000000_10 = std::string(1000000, '0');

    for (int i = 0; i < 1000; i += 10) {
        b_str_1000_10[i] = '1';
        if (i % 20 == 0) {
            b_str_1000_5[i] = '1';
            if (i % 100 == 0) {
                b_str_1000_1[i] = '1';
            }
        }
    }
    for (int i = 0; i < 10000; i += 10) {
        b_str_10000_10[i] = '1';
        if (i % 20 == 0) {
            b_str_10000_5[i] = '1';
            if (i % 100 == 0) {
                b_str_10000_1[i] = '1';
            }
        }
    }
    for (int i = 0; i < 100000; i += 10) {
        b_str_100000_10[i] = '1';
        if (i % 20 == 0) {
            b_str_100000_5[i] = '1';
            if (i % 100 == 0) {
                b_str_100000_1[i] = '1';
            }
        }
    }

    for (int i = 0; i < 1000000; i += 10) {
        b_str_1000000_10[i] = '1';
        if (i % 20 == 0) {
            b_str_1000000_5[i] = '1';
            if (i % 100 == 0) {
                b_str_1000000_1[i] = '1';
            }
        }
    }


    Bitvector b1000_1 = Bitvector(b_str_1000_1);
    Bitvector b1000_5 = Bitvector(b_str_1000_5);
    Bitvector b1000_10 = Bitvector(b_str_1000_10);
    Bitvector b10000_1 = Bitvector(b_str_10000_1);
    Bitvector b10000_5 = Bitvector(b_str_10000_5);
    Bitvector b10000_10 = Bitvector(b_str_10000_10);
    Bitvector b100000_1 = Bitvector(b_str_100000_1);
    Bitvector b100000_5 = Bitvector(b_str_100000_5);
    Bitvector b100000_10 = Bitvector(b_str_100000_10);

    Bitvector b1000000_1 = Bitvector(b_str_1000000_1);
    Bitvector b1000000_5 = Bitvector(b_str_1000000_5);
    Bitvector b1000000_10 = Bitvector(b_str_1000000_10);

    // Task 1
    // Jacobson's Rank
    // Timing
    // N = 1000, 10000, 100000

    Rank r = Rank(&b1000000_5);

    auto start = std::chrono::high_resolution_clock::now();

    r.rank1(10000);
    r.rank1(20000);
    r.rank1(30000);
    r.rank1(40000);
    r.rank1(50000);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    std::cout << r.overhead();

    // Task 2
    // lg(n)-time select operation that uses binary search over ranks

    Rank r = Rank(&b1000000_5);
    Select s = Select(&r);

    start = std::chrono::high_resolution_clock::now();

    s.select1(5);
    s.select1(10);
    s.select1(15);
    s.select1(20);
    s.select1(25);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    std::cout << s.overhead();

    // Task 3

    Sparse sp_1000_1 = Sparse(1000);
    Sparse sp_1000_5 = Sparse(1000);
    Sparse sp_1000_10 = Sparse(1000);
    //Sparse sp_10000_1 = Sparse(10000);
    //Sparse sp_10000_5 = Sparse(10000);
    //Sparse sp_10000_10 = Sparse(10000);
    //Sparse sp_100000_1 = Sparse(100000);
    //Sparse sp_100000_5 = Sparse(100000);
    //Sparse sp_100000_10 = Sparse(100000);
    //Sparse sp_1000000_1 = Sparse(1000000);
    //Sparse sp_1000000_5 = Sparse(1000000);
    //Sparse sp_1000000_10 = Sparse(1000000);

    for (int i = 0; i < 1000; i += 10) {
        sp_1000_10.append("foo",i);
        if (i % 20 == 0) {
            sp_1000_5.append("foo",i);
            if (i % 100 == 0) {
                sp_1000_1.append("foo",i);
            }
        }
    }
    /*for (int i = 0; i < 10000; i += 10) {
        //sp_10000_10.append("foo",i);
        if (i % 20 == 0) {
            sp_10000_5.append("foo",i);
            if (i % 100 == 0) {
                //sp_10000_1.append("foo",i);
            }
        }
    }
    for (int i = 0; i < 100000; i += 10) {
        //sp_100000_10.append("foo",i);
        if (i % 20 == 0) {
            sp_100000_5.append("foo",i);
            if (i % 100 == 0) {
                //sp_100000_1.append("foo",i);
            }
        }
    }

    for (int i = 0; i < 1000000; i += 10) {
        //sp_1000000_10.append("foo",i);
        if (i % 20 == 0) {
            sp_1000000_5.append("foo",i);
            if (i % 100 == 0) {
                //sp_1000000_1.append("foo",i);
            }
        }
    }*/

    sp_1000_1.finalize();
    sp_1000_5.finalize();
    sp_1000_10.finalize();
    //sp_10000_1.finalize();
    //sp_10000_5.finalize();
    //sp_10000_10.finalize();
    //sp_100000_1.finalize();
    //sp_100000_5.finalize();
    //sp_100000_10.finalize();
    //sp_1000000_1.finalize();
    //sp_1000000_5.finalize();
    //sp_1000000_10.finalize();

    std::string e = "hi";

    std::cout << "1000, 1" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000_1.get_at_rank(1, e);
    sp_1000_1.get_at_rank(2, e);
    sp_1000_1.get_at_rank(3, e);
    sp_1000_1.get_at_rank(4, e);
    sp_1000_1.get_at_rank(5, e);
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000_1.get_at_index(10, e);
    sp_1000_1.get_at_index(20, e);
    sp_1000_1.get_at_index(30, e);
    sp_1000_1.get_at_index(40, e);
    sp_1000_1.get_at_index(50, e);
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000_1.get_index_of(1);
    sp_1000_1.get_index_of(2);
    sp_1000_1.get_index_of(3);
    sp_1000_1.get_index_of(4);
    sp_1000_1.get_index_of(5);
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000_1.num_elem_at(10);
    sp_1000_1.num_elem_at(20);
    sp_1000_1.num_elem_at(30);
    sp_1000_1.num_elem_at(40);
    sp_1000_1.num_elem_at(50);
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;

/*    start = std::chrono::high_resolution_clock::now();

    sp_1000_10.get_at_rank(1, e);
    sp_1000_10.get_at_rank(2, e);
    sp_1000_10.get_at_rank(3, e);
    sp_1000_10.get_at_rank(4, e);
    sp_1000_10.get_at_rank(5, e);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "1000, 10" << std::endl;
    std::cout << duration.count() << std::endl;*/

/*    start = std::chrono::high_resolution_clock::now();

    sp_10000_1.get_at_rank(1, e);
    sp_10000_1.get_at_rank(2, e);
    sp_10000_1.get_at_rank(3, e);
    sp_10000_1.get_at_rank(4, e);
    sp_10000_1.get_at_rank(5, e);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "10000, 1" << std::endl;
    std::cout << duration.count() << std::endl;*/

/*    std::cout << "10000, 5" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_10000_5.size();
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_10000_5.num_elem();
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_10000_5.save("test.txt");
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_10000_5.load("test.txt");
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;*/

/*
    start = std::chrono::high_resolution_clock::now();

    sp_10000_10.get_at_rank(1, e);
    sp_10000_10.get_at_rank(2, e);
    sp_10000_10.get_at_rank(3, e);
    sp_10000_10.get_at_rank(4, e);
    sp_10000_10.get_at_rank(5, e);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "10000, 10" << std::endl;
    std::cout << duration.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();

    sp_100000_1.get_at_rank(1, e);
    sp_100000_1.get_at_rank(2, e);
    sp_100000_1.get_at_rank(3, e);
    sp_100000_1.get_at_rank(4, e);
    sp_100000_1.get_at_rank(5, e);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "100000, 1" << std::endl;
    std::cout << duration.count() << std::endl;
*/

/*    std::cout << "100000, 5" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_100000_5.size();
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_100000_5.num_elem();
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_100000_5.save("test.txt");
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_100000_5.load("test.txt");
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;*/

/*
    start = std::chrono::high_resolution_clock::now();

    sp_100000_10.get_at_rank(1, e);
    sp_100000_10.get_at_rank(2, e);
    sp_100000_10.get_at_rank(3, e);
    sp_100000_10.get_at_rank(4, e);
    sp_100000_10.get_at_rank(5, e);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "100000, 10" << std::endl;
    std::cout << duration.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();

    sp_1000000_1.get_at_rank(1, e);
    sp_1000000_1.get_at_rank(2, e);
    sp_1000000_1.get_at_rank(3, e);
    sp_1000000_1.get_at_rank(4, e);
    sp_1000000_1.get_at_rank(5, e);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "1000000, 1" << std::endl;
    std::cout << duration.count() << std::endl;
*/

/*    std::cout << "1000000, 5" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000000_5.size();
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000000_5.num_elem();
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000000_5.save("test.txt");
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sp_1000000_5.load("test.txt");
    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;*/

/*    start = std::chrono::high_resolution_clock::now();

    sp_1000000_10.get_at_rank(1, e);
    sp_1000000_10.get_at_rank(2, e);
    sp_1000000_10.get_at_rank(3, e);
    sp_1000000_10.get_at_rank(4, e);
    sp_1000000_10.get_at_rank(5, e);

    stop = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "1000000, 10" << std::endl;
    std::cout << duration.count() << std::endl;*/

    /*Sparse sp = Sparse(10);
    sp.append("foo", 1);
    sp.append("bar", 5);
    sp.append("baz", 9);
    sp.finalize();
    std::string e = "hi";
    //std::string e = "bar";
    std::cout << sp.get_at_rank(1, e) << std::endl;
    std::cout << e << std::endl;
    std::cout << sp.get_at_index(3, e) << std::endl;
    std::cout << e << std::endl;
    std::cout << sp.get_at_index(5, e) << std::endl;
    std::cout << e << std::endl;
    std::cout << sp.get_index_of(2);*/

    std::cout << sizeof (std::string) << std::endl;
    std::cout << "hi" << std::endl;

    return 0;
}

/* Scratch work
 */

/*
    for (const std::uint8_t i : { 0, 0b11111111, 0b00011101 }) {
        std::cout << "popcount( " << std::bitset<8>(i) << " ) = "
                  << std::popcount(i) << '\n';
        std::cout << sizeof(i) << std::endl;
    }*/

/*
 *     int popcounter = 0;
std::cout << std::bitset<100>(*ary.get()) << std::endl;
std::cout << std::bitset<64>(*ary.get()) << std::endl;
std::cout << std::bitset<64>(*(ary.get()+1)) << std::endl;
std::cout << std::popcount(*ary.get()) << std::endl;
std::cout << std::popcount(*(ary.get()+1)) << std::endl;
std::cout << popcounter << std::endl; */
/*
std::vector<int> ary1;
for (int k = 0; k < 25; k++) {
ary1.push_back(k);
}
std::cout << sizeof(std::vector<std::vector<int>>) << std::endl; */

/*
int my_func(int * t) {
    std::cout << *t << std::endl;
    return 1;
}*/

/*compact::vector<uint64_t, 1> ary(100);
std::string bit_vector_string;
//bit_vector_string = "0110101011011111010101010101101111000101111111341341341341341341341341341341341234";
//bit_vector_string = "101010101010101010101010101010101010101010101010101010101111111111111111111111111111111";
bit_vector_string = "100000101010101010101010101010101010101010101010101010101010101010101010101010101010101010101";
std::cout << bit_vector_string.size() << std::endl;
for (int i = 0; i < ary.size(); i++) {
    ary[i] = bit_vector_string[i];
    //std::cout << ary[i];
}
//std::cout << std::endl;
//std::cout << ary.size() << std::endl;
//std::cout << ary.bits() << std::endl;
std::cout << ary.bytes() << std::endl;
//std::cout << sizeof(uint64_t) << std::endl;



uint64_t * ptr = ary.get();
std:: cout << ptr << std::endl;
std::cout << ptr+1 << std::endl;
std::cout << ary.end() << std::endl;
int end = ary.bytes()/8;

for (int i = 0; i < end; i++) {
    std::cout << std::bitset<64>(*(ptr+i)) << std::endl;
}

Bitvector b = Bitvector(bit_vector_string);
std::cout << b.get_b_string() << std::endl;
std::cout << b.get_b_string().length() << std::endl;
std::cout << sizeof(b) << std::endl;
//int x = 5;
//my_func(&x);

//Rank r = Rank();
//std::cout << r.rank1(32) << std::endl;
//r.save("rank_save_test.txt");
//r.load("rank_save_test.txt");
//Select s = Select(&r);
//std::cout << r.rank1(128) << std::endl;
//std::cout << s.select1(45);
*/
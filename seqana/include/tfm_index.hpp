/*
 * tfm_index.hpp for BWT Tunneling
 * Copyright (c) 2020 Uwe Baier All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
//ORIGINALLY COMES FROM:
/*
 * tfm_index.hpp for Edge minimization in de Bruijn graphs
 * Copyright (c) 2019 Uwe Baier, Pascal Weber All Rights Reserved.
 */

#ifndef TFM_INDEX_HPP
#define TFM_INDEX_HPP

#include <sdsl/csa_wt.hpp>

#include <algorithm>
#include <limits>
#include <utility>
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;
using namespace sdsl;

//! a class representing a tunneled fm-index
template<class t_wt_type =       wt_blcd<>,
         class t_bv_type =       typename t_wt_type::bit_vector_type,
         class t_rank_type =     typename t_bv_type::rank_1_type,
         class t_select_type =   typename t_bv_type::select_1_type>
class tfm_index {
public:
	typedef int_vector<>::size_type           size_type;

	typedef typename t_wt_type::value_type          value_type;

	typedef t_wt_type                               wt_type;
    typedef t_bv_type                               bit_vector_type;
    typedef t_rank_type                             rank_type;
	typedef t_select_type                           select_type;

	//first index is next outgoing edge, second index is tunnel entry offset
	typedef pair<size_type,size_type>          nav_type;

private:
	template <typename t_tfm_index_type>
	friend void construct_tfm_index( t_tfm_index_type &tfm_index, uint64_t text_len, 
		int_vector_buffer<8> &&L_buf, bit_vector &&dout, bit_vector &&din );

	size_type                                       text_len; //original textlen
	wt_type                                         m_L;
	vector<size_type>                               m_C;
	bit_vector_type                                 m_dout;
	rank_type                                       m_dout_rank;
	select_type                                     m_dout_select;
    bit_vector_type                                 m_din;
	rank_type                                       m_din_rank;
	select_type                                     m_din_select;

public:
	const wt_type &                                        L = m_L;
	const vector<size_type> &                              C = m_C;
    const bit_vector_type &                                dout = m_dout;
	const rank_type &                                      dout_rank = m_dout_rank;
	const select_type &                                    dout_select = m_dout_select;
	const bit_vector_type &                                din = m_din;
	const rank_type &                                      din_rank = m_din_rank;
	const select_type &                                    din_select = m_din_select;

	//! returns the size of the original string
	size_type size() const {
		return text_len;
	};

	//! returns the end, i.e. the position in L where the string ends
	nav_type end() const {
		return make_pair( (size_type)0, (size_type)0 );
	}

	//! returns the character preceding the current position
/*	value_type preceding_char( const nav_type &pos ) const {
		return L[pos.first];
	}
*/
	//! Operation performs an backward step from current position.
	//! function sets pos to the new value and returns the result
	//! of preceding_char( pos ) before the backward step was performed
	pair <value_type, stack <size_type>> backwardstep( nav_type &pos, stack <size_type> s ) const {
		size_type &i = pos.first; //create references into position pair
		size_type &o = pos.second;

        //navigate to next entry
		auto is = L.inverse_select( i );
		auto c = is.second;
		i = C[c] + is.first;

		//check for the start of a tunnel
		auto din_rank_ip1 = din_rank( i + 1 );
		if (din[i] == 0 || din[i+1] == 0) {
			o = i - din_select( din_rank_ip1 ); //save offset to uppermost entry edge
            s.push(o);
		}

		//navigate to outedges of current node
		i = dout_select( din_rank_ip1 );

		//check for end of a tunnel
		if (dout[i+1] == 0) {
            o = s.top();
            i += o; //jump back offset
            s.pop();
		}

		return make_pair(c, s);
	};

	//! serializes opbject
	size_type serialize(ostream &out, structure_tree_node *v,
                          string name) const {

		structure_tree_node *child =
			structure_tree::add_child(v, name, util::class_name(*this));
		size_type written_bytes = 0;
		written_bytes += write_member(text_len, out, child, "text_len");

		written_bytes += m_L.serialize(out, child, "L");
		written_bytes += sdsl::serialize(m_C, out, child, "C");

		written_bytes += m_dout.serialize(out, child, "dout");
		written_bytes += m_dout_rank.serialize(out, child, "dout_rank");
		written_bytes += m_dout_select.serialize(out, child, "dout_select");

		written_bytes += m_din.serialize(out, child, "din");
		written_bytes += m_din_rank.serialize(out, child, "din_rank");
		written_bytes += m_din_select.serialize(out, child, "din_select");

		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	};

	//! loads a serialized object
	void load(istream &in) {

		read_member( text_len, in );

		m_L.load(in);
		sdsl::load(m_C, in);

		m_dout.load(in);
		m_dout_rank.load(in, &m_dout);
		m_dout_select.load(in, &m_dout);

		m_din.load(in);
		m_din_rank.load(in, &m_din);
		m_din_select.load(in, &m_din);
	};
};

//// SPECIAL CONSTRUCTION FOR TUNNELED FM INDEX ///////////////////////////////

//! Function computes, whether there is an element in the column that is also in the presorted block using binary search
//! Function return true if there exists such an element
template <class t_tfm_index_type>
bool column_in_sorted_block (vector <typename t_tfm_index_type::size_type> column,
                               vector <typename t_tfm_index_type::size_type> sorted_block, typename t_tfm_index_type::size_type width,
                               typename t_tfm_index_type::size_type shift) {
    typedef typename t_tfm_index_type::size_type size_type;

    for (size_type i = 0; i < column.size(); i++) {
        size_type left = 0;
        size_type right = sorted_block.size();

        while (left < right) {
            size_type tmp = (left + right)/2;

            if (column[i] + shift >= sorted_block[tmp] && column[i] + shift < sorted_block[tmp] + width)
                return true;
            else if (column[i] + shift >= sorted_block[tmp])
                left = tmp;
            else
                right = tmp;

            if (left >= right-1 && left == tmp)
                break;
        }
    }

    return false;
}

//! Function computes, whether there is an element in the presorted column that is also in the block using binary search
//! Function return true if there exists such an element
template <class t_tfm_index_type>
bool sorted_column_in_block (vector <typename t_tfm_index_type::size_type> column, vector <typename t_tfm_index_type::size_type> block,
                               typename t_tfm_index_type::size_type width, typename t_tfm_index_type::size_type shift) {
    typedef typename t_tfm_index_type::size_type size_type;

    for (size_type i = 0; i < block.size(); i++) {
        size_type left = 0;
        size_type right = column.size();

        while (left < right) {
            size_type tmp = (left + right)/2;

            if ((column[tmp] + shift) >= block[i] && (column[tmp] + shift) < (block[i] + width))
                return true;
            else if ((column[tmp] + shift) < block[i])
                left = tmp;
            else
                right = tmp;

            if (left >= right-1 && left == tmp)
                break;
        }
    }
    return false;
}

//! Function enumerates all height-maximal, left-maximal non-self-colliding blocks in the csa
//! and computes all collisions between them, storing them in vectors critical_collisions and compensable_collisions
//! Function returns vector blocks with blocks stored as {w, {a, b}},
//! where w is width of the block and <a, b) is half-closed interval of the first column of the block
template <class t_tfm_index_type, class t_csa_wt_type>
vector <pair <typename t_tfm_index_type::size_type, pair <typename t_tfm_index_type::size_type, typename t_tfm_index_type::size_type>>>
        find_blocks_and_collisions (t_csa_wt_type csa,
                                    vector <pair <typename t_tfm_index_type::size_type, typename t_tfm_index_type::size_type>> &critical_collisions,
                                    vector <pair <typename t_tfm_index_type::size_type, typename t_tfm_index_type::size_type>> &compensable_collisions) {

    typedef typename t_tfm_index_type::size_type size_type;

    vector <pair <size_type, pair <size_type, size_type>>> all_blocks;

    //compute LCS array
    int_vector<> LCS(csa.size(), 0);
    {
        vector<size_type> LF_inverse (csa.lf.size());
        for (size_type i = 0; i < csa.lf.size(); i++)
            LF_inverse[csa.lf[i]] = i;

        LCS[1] = 0;
        size_type j = LF_inverse[1];
        size_type l = 0;
        for (size_type i = 2; i < csa.bwt.size(); i++) {
            if (j > 0 && csa.bwt[j] == csa.bwt[j-1])
                l++;
            else
                l = 0;
            LCS[j] = l;
            j = LF_inverse[j];
        }
    }

    //find all left-maximal and height-maximal blocks
    // store beginnings of potential blocks as <start position, width> on stack
    stack <pair <size_type, size_type>> stack;
    stack.push(make_pair(1, 0));

    for (size_type i = 1; i < LCS.size(); i++) {
        pair <size_type, size_type> start (stack.top());

        //end of the block of width of size start.second
        while (start.second > LCS[i]) {
            stack.pop();

            //if big enough, report the block
            if (start.second > 1 && i - start.first > 1 /*&& start.second%2 == 0*/) {
                all_blocks.push_back (make_pair (start.second+1, make_pair (csa.isa[csa[start.first]-start.second], csa.isa[csa[i-1]-start.second]+1)));
            }

            //assurance of the height-maximality of blocks
            if (LCS[i] > 1 && stack.top().second < LCS[i])
                stack.push (make_pair (start.first, LCS[i]));
            start = stack.top();
        }

        //store a possible start of a block
        if (start.second < LCS[i])
            stack.push (make_pair (i-1, LCS[i]));
    }

    while (!stack.empty()) {
        pair <size_type, size_type> start (stack.top());
        if (start.second > 1 && csa.size()-start.first > 1)
            all_blocks.push_back (make_pair (start.second+1, make_pair (csa.isa[csa[start.first]-start.second], csa.isa[csa[csa.size()-1]-start.second]+1)));
        stack.pop();
    }

    //get rid of self-colliding blocks
    vector <pair <size_type, pair <size_type, size_type>>> blocks;
    for (pair <size_type, pair <size_type, size_type>> block : all_blocks) {
        bool not_self_colliding = true;

        vector <size_type> start_pos(csa.begin() + block.second.first, csa.begin() + block.second.second);
        sort(start_pos.begin(), start_pos.end());

        for (size_type m = 1; m < start_pos.size(); m++) {
            if (start_pos[m] - start_pos[m-1] < block.first) {
                not_self_colliding = false;
                break;
            }
        }

        //store the block if it is non-self-colliding and determine th type of collision with previously checked blocks
        if (not_self_colliding) {
            for (size_type i = 0; i < blocks.size(); i++) {
                pair <size_type, pair <size_type, size_type>> block2 = blocks[i];
                vector <size_type> start_pos2 (csa.begin() + block2.second.first, csa.begin() + block2.second.second);

                for (int i = 0; i < 2; i++) {
                //check if it is an obvious critical collision (first positions is equal or first column of the wider block is shared)
                if ((block2.second.first == block.second.first) || (block2.first >= block.first && column_in_sorted_block<t_tfm_index_type> (start_pos2, start_pos, block.first, 0))
                    || (block2.first <= block.first && sorted_column_in_block<t_tfm_index_type> (start_pos, start_pos2, block2.first, 0))) {
                    critical_collisions.push_back(make_pair(i, blocks.size()));
                    continue;
                }
                }
                /*
                //check if they do collide, i.e. if the first positions of the inner block are not shared
                if ((block2.first > block.first && !sorted_column_in_block<t_tfm_index_type> (start_pos, start_pos2, block2.first, 0))
                    || (block2.first <= block.first && !column_in_sorted_block<t_tfm_index_type> (start_pos2, start_pos, block.first, 0)))
                    continue;

                //indicate how blocks collide (if the last column of wider block is not shared, it is compensable collision)
                if ((block2.first < block.first && !sorted_column_in_block<t_tfm_index_type> (start_pos, start_pos2, block2.first, block.first - 1))
                    || (block2.first > block.first && !column_in_sorted_block<t_tfm_index_type> (start_pos2, start_pos, block.first, block2.first - 1))) {
                    if (block2.first > block.first) {
                        compensable_collisions.push_back(make_pair(blocks.size(), i));
                    } else {
                        compensable_collisions.push_back(make_pair(i, blocks.size()));
                    }
                } else {
                    critical_collisions.push_back(make_pair(i, blocks.size()));
                }*/
            }
            blocks.push_back(block);
        }
    }
    return blocks;
}


//! Function constructs a tfm index using Gurobi ilp solver
//! Function returns the price of the compression
template <class t_tfm_index_type, class t_csa_wt_type>
string construct_tbwt ( t_tfm_index_type &tfm_index, t_csa_wt_type &&csa, cache_config &config ) {
    typedef typename t_tfm_index_type::size_type size_type;
    string result;

    //find all non self-colliding blocks and store them in vector blocks as <width, interval in sa>,
    // store critical collisions between blocks in vector critical_collisions
    // and compensable collisions in vector compensable_collisions
    vector <pair<size_type, size_type>> critical_collisions, compensable_collisions;
    vector <pair <size_type, pair <size_type, size_type>>> blocks =
            find_blocks_and_collisions <t_tfm_index_type, t_csa_wt_type> (csa, critical_collisions, compensable_collisions);

    //compute block price and price of shared block in compensable collision
    vector<long double> block_price, shared_price;
    {
        for (pair <size_type, pair<size_type, size_type>> block: blocks) {
            long double price = (block.first - 1)* (block.second.second - block.second.first - 1);
            block_price.push_back(price);
        }

        for (pair <size_type, size_type> cc: compensable_collisions) {
            long double price = (blocks[cc.first].first - 1) * (blocks[cc.second].second.second - blocks[cc.second].second.first - 1);
            shared_price.push_back(price);
        }
    }

    //creating ilp instance
    {
        ofstream outfile("reduction.lp");

        //goal and main function
        outfile << "Maximize" << endl;
        for (size_type i = 0; i < blocks.size(); i++) {
            if (i > 0)
                outfile << " + ";
            outfile << block_price[i] << " x" << i;
        }
        for (size_type i = 0; i < shared_price.size(); i++) {
            outfile << " - " << shared_price[i] << " y" << i;
        }

        //constraints
        outfile << endl << "Subject To" << endl;

        //constraints for critical collisions
        for (size_type i = 0; i < critical_collisions.size(); i++) {
            outfile << "crc" << i << ": x" << critical_collisions[i].first << " + x" << critical_collisions[i].second
                    << " <= 1" << endl;
        }

        //constraints for compensable collisions
        vector <vector <size_type>> comp_coll (blocks.size());
        for (pair <size_type, size_type> coc : compensable_collisions)
            comp_coll[coc.first].push_back(coc.second);

        for (size_type i = 0; i < compensable_collisions.size(); i++) {
            size_type Bin = compensable_collisions[i].first;
            size_type Bout = compensable_collisions[i].second;

            vector<size_type> inter_blocks;

            //remember compensable collisions as vector of outer neighbours
            for (size_type j : comp_coll[Bin]) {
                if (j != Bout) {
                    for (size_type k : comp_coll[j]) {
                        if (k == Bout)
                            inter_blocks.push_back(j);
                    }
                }
            }

            outfile << "coc" << i << "-0: x" << Bin<< " + x"
                    << Bout << " + " << inter_blocks.size();
            for (size_type ib : inter_blocks)
                outfile << " - x" << ib;
            outfile << " - " << inter_blocks.size()+2 << " y" << i << " >= 0" << endl;

            outfile << "coc" << i << "-1: x" << Bin << " + x" << Bout;
            for (size_type ib : inter_blocks)
                outfile << " - x" << ib;
            outfile << " - " << inter_blocks.size()+2 << " y" << i << " <= 1" << endl;
        }

        //boundaries
        outfile << "Binaries" << endl;
        outfile << "x0";
        for (size_type i = 1; i < blocks.size(); i++) {
            outfile << " x" << i;
        }
        for (size_type i = 0; i < compensable_collisions.size(); i++) {
            outfile << " y" << i;
        }

        outfile << endl << "End" << endl;
        outfile.close();
    }

    system("./../gurobi9.5.1_linux64/gurobi951/linux64/bin/gurobi_cl ResultFile=reduction.sol reduction.lp >gurobi_output.txt");
    cout << "The output from gurobi solver can be found in gurobi_output.txt" << endl;

    //load blocks for tunneling from the Gurobi solution file
    bit_vector tunneled(blocks.size(), 0);
    {
        ifstream solution;
        solution.open("reduction.sol");

        if (!solution) { // file couldn't be opened
            cerr << "Error: file reduction.sol could not be opened" << endl;
            exit(1);
        }

        string str;
        solution >> str >> str >> str >> str >> result;
        for (size_type i = 0; !solution.eof() && i < blocks.size(); i++) {
            int tunnel;
            solution >> str >> tunnel;
            if (tunnel == 1)
                tunneled[i] = 1;
        }
        solution.close();
    }

    //tunneling bwt

    //marking to-be-tunneled blocks in arrays din and dout
    {
        bit_vector dout (csa.bwt.size() + 1, 1);
        bit_vector din (csa.bwt.size() + 1, 1);
        for (size_type i = 0; i < tunneled.size(); i++) {
            if (tunneled[i] != 1)
                continue;

            size_type a = blocks[i].second.first;
            size_type h = blocks[i].second.second - a;

            for (size_type j = 0; j < blocks[i].first - 1; j++) {
                for (size_type k = 1; k < h; k++)
                    dout[a + k] = 0;

                a = csa.isa[csa[a] + 1];

                for (size_type k = 1; k < h; k++)
                    din[a + k] = 0;
            }
        }

        //create a buffer for newly constructed L
        string tmp_key = util::to_string (util::pid()) + "_" + util::to_string (util::id());
        string tmp_file_name = cache_file_name (tmp_key, config);
        int_vector_buffer<8> L_buf (tmp_file_name, ios::out);

        //remove redundant entries from L, dout and din
        size_type p = 0, q = 0;
        for (size_type i = 0; i <  csa.bwt.size(); i++) {
            if (din[i] == 1) {
                L_buf.push_back(csa.bwt[i]);
                dout[p++] = dout[i];
            }

            if (dout[i] == 1) {
                din[q++] = din[i];
            }
        }
        dout[p++] = 1;
        din[q++] = 1;
        dout.resize(p);
        din.resize(q);

        uint64_t text_len = csa.size();
        csa = t_csa_wt_type(); //remove csa object as it is no longer required
	
        cout << "L" << endl;
        for (size_type i = 0; i < L_buf.size(); i++)
                cout << L_buf[i];
        cout << endl << endl;


	cout << "dout" << endl;
	for (size_type i = 0; i < p; i++)
		cout << dout[i];
	cout << endl << endl;
	
	cout << "din" << endl;
	for (size_type i = 0; i < q; i++)
		cout << din[i];
	cout << endl;

        construct_tfm_index(tfm_index, text_len, move(L_buf), move(dout), move(din));

        remove(tmp_file_name);
    }

    return result;
};

//! Function constructs tfm_index from L, din and dout
template<class t_tfm_index_type>
void construct_tfm_index( t_tfm_index_type &tfm_index, uint64_t text_len, int_vector_buffer<8> &&L_buf, bit_vector &&dout, bit_vector &&din ) {
    //set original string size
    tfm_index.text_len = text_len;

    //construct tfm index from L, din and dout
    typedef typename t_tfm_index_type::wt_type         wt_type;
    typedef typename t_tfm_index_type::bit_vector_type bv_type;

    //wavelet tree of L
    tfm_index.m_L = wt_type( L_buf, L_buf.size() );
    create_C_array( tfm_index.m_C, tfm_index.m_L );

    //dout
    tfm_index.m_dout = bv_type( move( dout ) );
    util::init_support( tfm_index.m_dout_rank, &tfm_index.m_dout );
    util::init_support( tfm_index.m_dout_select, &tfm_index.m_dout );

    //din
    tfm_index.m_din = bv_type( move( din ) );
    util::init_support( tfm_index.m_din_rank, &tfm_index.m_din );
    util::init_support( tfm_index.m_din_select, &tfm_index.m_din );
};

#endif

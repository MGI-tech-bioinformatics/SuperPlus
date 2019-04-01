// MakeDepend: library HTSLIB
// MakeDepend: cflags HTSLIB_CFLAGS
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include "CoreTools.h"
#include "MainTools.h"
#include "FastIfstream.h"
#include "ParallelVecUtilities.h"
#include "10X/Barcode.h"
#include "system/System.h"
#include <htslib/sam.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <algorithm>
#include <ctime>

namespace dr_wu {
	void parse_paf(const char* paf_file, const vec<String>& scaff_name, vec<Bool>& non_dup)
	{
		std::cerr << Date() << ": start to parse paf file" << std::endl;
		non_dup.clear();
		non_dup.resize(scaff_name.size(), true);


		std::ifstream paf_stream;
		paf_stream.open(paf_file);
		if (paf_stream.fail())
		{
			std::cerr << Date() << ": Can't read paf file, continue." << std::endl;
			return;
		}
		std::string scaff1, scaff2;
		int length1, length2;
		char dir;
		int qstart, qend;
		int tstart, tend;
		int rmatch;
		int align_length;
		int mapq;
		std::string sam;
		std::string line;
		std::stringstream linestream;
		vec<std::pair<std::string, std::string>> overlaps;
		std::map<std::string, int> lengths_map;
		std::set<std::string> scaffs;

		std::set<std::string> non_dup_scaffs;

		while (!paf_stream.eof())
		{
			std::getline(paf_stream, line);
			linestream.str(line);
			linestream >> scaff1 >> length1 >> qstart >> qend >> dir >> scaff2 >> length2 >> tstart >> tend >> rmatch >> align_length >> mapq >> sam;
			lengths_map[scaff1] = length1;
			lengths_map[scaff2] = length2;
			scaffs.insert(scaff1);
			scaffs.insert(scaff2);
			if (scaff1 == scaff2) continue;
			if (align_length > length1 * 9 / 10) { overlaps.push_back(std::make_pair(scaff1, scaff2)); }
			if (align_length > length2 * 9 / 10) { overlaps.push_back(std::make_pair(scaff2, scaff1)); }
			//std::cerr<<scaff1<<"\t"<<scaff2<<"\t"<<sam<<std::endl;
		}

		paf_stream.close();

		std::map<std::string, int> scaff_n;
		vec<std::string> scaff_v;
		vec<int> scaff_length;
		int i = 1;
		for (auto& scaff : scaffs)
		{
			scaff_v.push_back(scaff);
			scaff_n[scaff] = i;
			scaff_length.push_back(lengths_map[scaff]);
			i++;
		}

		int count = scaffs.size();
		vec<std::set<int>> to_(count);
		vec<std::set<int>> from_(count);

		for (auto& scaff_pair : overlaps)
		{
			int i1 = scaff_n[scaff_pair.first] - 1;
			int i2 = scaff_n[scaff_pair.second] - 1;
			from_[i1].insert(i2);
			to_[i2].insert(i1);
		}

		vec<Bool> proceeded(count, false); //mark elements not maximal
		vec<Bool> maximal(count, false); //mark maximal elements

		for (int i = 0; i < count; i++)
		{
			if (from_[i].size() == 0)
			{
				maximal[i] = true;
				//std::cerr<<scaff_v[i]<<"\t"<<scaff_length[i]<<std::endl;
				non_dup_scaffs.insert(scaff_v[i]);
				continue;
			}
			if (maximal[i] || proceeded[i]) { continue; }
			std::set<int> s1 = {}; //set of smaller elements
			std::set<int> s2 = {}; //set of greater elements 

			//cerr<<from_[i].size()<<to_[i].size()<<endl;

			//use DFS to obtain s1&s2
			std::stack<int> stack1;
			stack1.push(i);
			while (!stack1.empty())
			{
				int v1 = stack1.top();
				stack1.pop();
				if (s1.find(v1) == s1.end())
				{
					s1.insert(v1);
					for (auto it = to_[v1].begin(); it != to_[v1].end(); it++)
					{
						stack1.push(*it);
					}
				}
			}
			stack1.push(i);
			while (!stack1.empty())
			{
				int v1 = stack1.top();
				stack1.pop();
				if (s2.find(v1) == s2.end())
				{
					s2.insert(v1);
					for (auto it = from_[v1].begin(); it != from_[v1].end(); it++)
					{
						stack1.push(*it);
					}
				}
			}

			/*cerr<<scaff_list[i]<<endl;
			  for(auto it=s1.begin();it!=s1.end();it++)
			  {
			  cerr<<scaff_list[*it]<<"\t";
			  }
			  cerr<<endl;
			  cerr<<endl;*/

			  //if s2 \subset s1, all elements in s2 are maximal, s1\setminus s2 not.
			  //else all element in s1 not.
			Bool issubset = true;
			for (auto it = s2.begin(); it != s2.end(); it++)
			{
				if (s1.find(*it) == s1.end()) { issubset = false; break; }
			}
			if (issubset)
			{
				int max_length = 0;
				int max_scaff = 0;
				for (auto it = s2.begin(); it != s2.end(); it++)
				{
					maximal[*it] = true;
					if (max_length < scaff_length[*it]) { max_length = scaff_length[*it]; max_scaff = *it; }
				}
				//std::cerr<<scaff_v[i]<<"\t"<<scaff_length[i]<<std::endl;

				non_dup_scaffs.insert(scaff_v[i]);
				//cerr<<scaff_list[max_scaff]<<endl;
				for (auto it = s1.begin(); it != s1.end(); it++)
				{
					if (s2.find(*it) == s2.end()) { proceeded[*it] = true; }
				}
			}
			else
			{
				for (auto it = s1.begin(); it != s1.end(); it++)
				{
					proceeded[*it] = true;
				}
			}
		}
		paf_stream.close();
		//for(int i=0;i<overlaps.size();i++){std::cerr<<overlaps[i].first<<"\t"<<overlaps[i].second<<std::endl;}

		for (int i = 0; i < scaff_name.isize(); i++)
		{
			non_dup[i] = (scaffs.find(scaff_name[i]) == scaffs.end()) || (non_dup_scaffs.find(scaff_name[i]) != scaffs.end());
		}
		std::cerr << Date() << ": end of parsing paf file" << std::endl;
	}

	int translate_coord(int coord, const vec<std::pair<int, int>>& largeN)
	{
		if (largeN.size() == 0) return coord;
		else
		{
			int acc_length = largeN[0].first, loop = 0;
			while (loop + 1 < largeN.isize() && acc_length < coord)
			{
				acc_length += largeN[loop + 1].first - largeN[loop].second;
				loop++;
			}
			//cerr<<"loop:"<<loop<<endl;
			if (acc_length > coord)
			{
				//cerr<<"case 1:"<<largeN[loop].first<<"\t"<<(acc_length-coord)<<endl;
				return largeN[loop].first - (acc_length - coord);
			}
			else
			{
				//cerr<<"case 2:"<<largeN[loop].second<<"\t"<<(acc_length-coord)<<endl;
				return largeN[loop].second - (acc_length - coord);
			}
		}
	}

	template <class T>
	int intersection_size(const std::set<T>& a, const std::set<T>& b)
	{
		int i = 0;
		for (auto s : a) { if (b.find(s) != b.end()) i++; }
		return i;
	}

	vec<vec<int> > star(const vec<quad<int, int, int, double> >& jac)
	{
		vec<vec<int>> new_scaffs;
		Bool verbose = false;

		std::map<std::pair<int, int>, double> adj_matrix;
		std::map<int, vec<std::pair<int, double> > > adj_list;
		std::map<int, double> max_jac;
		for (auto it = jac.begin(); it != jac.end(); it++)
		{
			adj_matrix[std::make_pair(it->first, it->second)] = it->fourth;
			adj_matrix[std::make_pair(it->second, it->first)] = it->fourth;
			if (max_jac[it->first] < it->fourth) { max_jac[it->first] = it->fourth; }
			if (max_jac[it->second] < it->fourth) { max_jac[it->second] = it->fourth; }
		}

		for (auto it = adj_matrix.begin(); it != adj_matrix.end(); it++)
		{
			adj_list[it->first.first].push_back(std::make_pair(it->first.second, it->second));
		}


		std::map<int, Bool> visited;
		vec<std::pair<int, int>> stars;
		std::map<int, int> linear_link;

		//this is only picking linear part!

		if (verbose) std::cerr << "linear links:" << std::endl;
		for (auto it = adj_list.begin(); it != adj_list.end(); it++)
		{
			if (it->second.size() == 1)
			{
				int s1 = it->second[0].first;
				if (adj_list[s1].size() == 1 && it->first > s1)
				{
					if (verbose) std::cerr << it->first << "\t" << it->second[0].first << std::endl;
					linear_link[s1] = it->first;
					linear_link[it->first] = s1;
				}
			}
		}

		std::cerr << "make " << linear_link.size() << " links" << std::endl;
		//cerr<<"scaffolds:"<<endl<<endl;
		int i = 0;
		for (auto it = linear_link.begin(); it != linear_link.end(); it++)
		{
			if (!visited[it->first / 2] && it->second != 0 && linear_link[it->first ^ 1] == 0)
			{

				new_scaffs.push_back({});
				int s = it->first ^ 1;
				visited[s / 2] = true;
				//cerr<<"Begin point:"<<s <<"\t next:"<<linear_link[s^1]<<endl;
				while (linear_link[s ^ 1] != 0)
				{
					//cerr<<s/2<<"\t"<<s%2<<endl;
					new_scaffs[i].push_back(s);
					s = linear_link[s ^ 1];
					visited[s / 2] = true;
				}

				new_scaffs[i].push_back(s);
				//cerr<<s/2<<"\t"<<s%2<<endl;
				//cerr<<endl;
				i++;
			}
		}


		return new_scaffs;
	}

	Bool sort_by_sec(const std::pair<int, int> a, const std::pair<int, int> b) { return a.second > b.second; }
	Bool sort_by_first(const std::pair<int, int> a, const std::pair<int, int> b) { return a.first > b.first; }

	int max(int a, int b) { return a > b ? a : b; }

	void linear_scaffolding(const vec<vec<std::pair<int, int> > >& lbp,
		const vec<int>& lengths,
		const vec<vec<std::pair<int, int> > >& largeNs,
		const vec<Bool>& non_dup,
		const vec<String>& scaff_list,
		vec<vec<int> >& new_scaffs, int window)
	{
		std::cerr << Date() << ": start to rescaffold" << std::endl;
		int FLANKING_REGION = window;
		const int MIN_READS = 5;
		const double MIN_JAC = 0.01;
		const int MIN_COMMON_BC = 20;
		const double MIN_LONG_FOLD = 10;
		const int MIN_LONG_BC = 10;

		//for()

		vec<Bool> valid(lengths.size(), true);
		vec<int> lengths_noN(lengths.size());

		vec<std::set<int>> scaff_bc_all(lengths.size());
		vec<std::set<int>> scaff_bc(lengths.size() * 2);
		vec<std::set<int>> scaff_bc_long(lengths.size() * 2);

#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < lengths.isize(); i++)
		{
			if (!non_dup[i]) { valid[i] = false; continue; }
			if (lengths[i] < 2 * FLANKING_REGION) { valid[i] = false; continue; }
			int length = lengths[i];
			for (int j = 0; j < largeNs[i].isize(); j++) { length -= (largeNs[i][j].second - largeNs[i][j].first); }
			if (length < 2 * FLANKING_REGION) { valid[i] = false; continue; }
			lengths_noN[i] = length;
			auto& largeN = largeNs[i];

			int l1 = translate_coord(FLANKING_REGION, largeN);
			int l2 = translate_coord(FLANKING_REGION * 2, largeN);
			int l3 = translate_coord(lengths_noN[i] - FLANKING_REGION, largeN);
			int l4 = translate_coord(lengths_noN[i] - FLANKING_REGION * 2, largeN);

			auto& bps = lbp[i];

			std::map<int, int> r_count_b;
			std::map<int, int> r_count_e;
			std::map<int, int> r_count_bl;
			std::map<int, int> r_count_el;
			std::map<int, int> r_count_total;

			for (auto& bp : bps)
			{
				r_count_total[bp.second]++;
				if (bp.first < l1) { r_count_b[bp.second]++; }
				else if (bp.first < l2) { r_count_bl[bp.second]++; }
				if (bp.first > l3) { r_count_e[bp.second]++; }
				else if (bp.first > l4) { r_count_el[bp.second]++; }
			}

			for (auto& bc : r_count_b) { if (bc.second > MIN_READS) { scaff_bc[2 * i].insert(bc.first); } }
			for (auto& bc : r_count_e) { if (bc.second > MIN_READS) { scaff_bc[2 * i + 1].insert(bc.first); } }
			for (auto& bc : r_count_bl) { if (bc.second > MIN_READS) { scaff_bc_long[2 * i].insert(bc.first); } }
			for (auto& bc : r_count_el) { if (bc.second > MIN_READS) { scaff_bc_long[2 * i + 1].insert(bc.first); } }
			for (auto& bc : r_count_total) { if (bc.second > MIN_READS) { scaff_bc_all[i].insert(bc.first); } }
		}

		vec<quad<int, int, int, double> > total_matrix;
		quad<int, int, int, double>  tmp;

		for (int i = 0; i < lengths.isize() * 2; i++)
		{
			if (!valid[i / 2]) continue;
			for (int j = 0; j < i; j++)
			{
				if (!valid[j / 2]) continue;
				auto& bcs1 = scaff_bc[i];
				auto& bcs2 = scaff_bc[j];
				auto& bcs_long1 = scaff_bc_long[i];
				auto& bcs_long2 = scaff_bc_long[j];
				if (i / 2 == j / 2)continue;
				if (intersection_size(bcs1, bcs2) < MIN_COMMON_BC) continue;
				double jac = (double)intersection_size(bcs1, bcs2) / (bcs1.size() + bcs2.size() - intersection_size(bcs1, bcs2));
				if (jac < MIN_JAC) continue;
				//cerr<<intersection_size(bcs_long1,bcs_long2)<<endl;
				if (intersection_size(bcs_long1, bcs_long2) < MIN_LONG_BC && MIN_LONG_FOLD*intersection_size(bcs_long1, bcs_long2) < intersection_size(bcs1, bcs2)) continue;

				//cal span
				std::set<int> common_all;
				for (auto& b : scaff_bc_all[i / 2]) { if (scaff_bc_all[j / 2].find(b) != scaff_bc_all[j / 2].end()) common_all.insert(b); }

				vec<int> pos1, pos2;

				/*for(auto& b:common_all)
				  {
				  pos1.insert(pos1.end(),lbp[i/2][b].begin(),lbp[i/2][b].end());
				  pos2.insert(pos2.end(),lbp[j/2][b].begin(),lbp[j/2][b].end());
				  }*/

				for (auto bp : lbp[i / 2])
				{
					if (common_all.find(bp.second) != common_all.end()) { pos1.push_back(bp.first); }
				}


				for (auto bp : lbp[j / 2])
				{
					if (common_all.find(bp.second) != common_all.end()) { pos2.push_back(bp.first); }
				}

				sort(pos1.begin(), pos1.end());
				sort(pos2.begin(), pos2.end());

				int m_pos1 = pos1[pos1.size() / 2];
				int m_pos2 = pos2[pos2.size() / 2];


				// std::cerr<<lengths[i/2]<<"\t"<<lengths[j/2]<<"\t"<<m_pos1<<"\t"<<m_pos2<<std::endl;

				if (i % 2 == 1) { m_pos1 = lengths[i / 2] - m_pos1; }

				if (j % 2 == 1) { m_pos2 = lengths[j / 2] - m_pos2; }

				//cerr<<length_v[i/2]<<"\t"<<length_v[j/2]<<"\t"<<m_pos1<<"\t"<<m_pos2<<endl;

				//if(m_pos2>1000000 || m_pos1>1000000) continue;

				tmp.first = i;
				tmp.second = j;
				tmp.third = intersection_size(bcs1, bcs2);
				tmp.fourth = jac;
#pragma omp critical(total_matrix)
				{
					total_matrix.push_back(tmp);
				}
				std::cerr << scaff_list[i / 2] << "\t" << i % 2 + 1 << "\t" << scaff_list[j / 2] << "\t" << j % 2 + 1 << "\t" << intersection_size(bcs1, bcs2) << "\t" << jac << std::endl;
			}
		}

		for (int i = 0; i < total_matrix.isize(); i++)
		{
			for (int j = 0; j < total_matrix.isize(); j++)
			{
				if (i == j) continue;
				if (total_matrix[i].first / 2 != total_matrix[j].first / 2 || total_matrix[i].second / 2 != total_matrix[j].second / 2) continue;
				if (total_matrix[i].fourth / 2 > total_matrix[j].fourth / 2) { total_matrix.erase(total_matrix.begin() + j); j--; }
			}
		}

		new_scaffs = star(total_matrix);

		std::set<int> linked;
		for (auto& new_scaff : new_scaffs)
		{
			for (auto& scaff : new_scaff)
			{
				linked.insert(scaff / 2);
				std::cout << scaff_list[scaff / 2] << "\t" << scaff % 2 << std::endl;
			}
			std::cout << std::endl;
		}

		for (int i = 0; i < lengths.isize(); i++)
		{
			if (linked.find(i) == linked.end() && non_dup[i])
			{
				new_scaffs.push_back({ 2 * i });
			}
		}
		std::cerr << Date() << ": end of rescaffolding" << std::endl;
	}

}

void ParseBAM(const char *bam_input,
	vec<vec<pair<int, int> > > &lbps,
	vec<int> &scaffold_lens,
	vec<String> &scaffold_names) {

	std::cerr << Date() << ": start to parse bam file" << std::endl;
	samFile *bam_fp = nullptr;
	bam1_t *bam_record = bam_init1();

	// open bam file.
	if ((bam_fp = sam_open(bam_input, "r")) == nullptr) {
		std::cerr << Date() << ": Can't open bam file: " << bam_input << std::endl;
		Scram(1);
	}

	// get header info.
	bam_hdr_t *bam_header = nullptr;
	if ((bam_header = sam_hdr_read(bam_fp)) == nullptr) {
		std::cerr << Date() << ": Can't read bam header." << std::endl;
		Scram(1);
	}

	// clear lbp
	lbps.clear_and_resize(bam_header->n_targets);

	scaffold_lens.clear();
	scaffold_names.clear();
	// record scaffold length and name
	for (int i = 0; i < bam_header->n_targets; ++i) {
		scaffold_lens.push_back(bam_header->target_len[i]);
		scaffold_names.push_back(bam_header->target_name[i]);
	}

	int64_t read_cnt = 0, read_filtered = 0;
	while (sam_read1(bam_fp, bam_header, bam_record) >= 0) {
		++read_cnt;
		if (bam_record->core.qual < 60) {
			++read_filtered;
			continue;
		}
		if (bam_record->core.flag & (0x4 | 0x100 | 0x200 | 0x400 | 0x800)) {
			++read_filtered;
			continue;
		}
		String barcode_str = reinterpret_cast<char*>(bam_record->data);
		barcode_str = barcode_str.SafeAfterLast("#");
		barcode_t barcode_int = createBarcodeIndex(barcode_str);
		if (barcode_int == 0) {
			++read_filtered;
			continue;
		}
		// save barcode locations to lbp
		lbps[bam_record->core.tid].push(
			((bam_endpos(bam_record) + bam_record->core.pos) / 2),
			barcode_int);
	}

	// sort lbp
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < lbps.isize(); ++i) {
		Sort(lbps[i]);
	}

	std::cerr << Date() << ": total reads: " << read_cnt << std::endl;
	std::cerr << Date() << ": filtered reads: " << read_filtered << std::endl;
	std::cerr << Date() << ": qualified reads: " << read_cnt - read_filtered
		<< std::endl;
	std::cerr << Date() << ": end of parsing bam file" << std::endl;
	bam_hdr_destroy(bam_header);
	bam_destroy1(bam_record);
	sam_close(bam_fp);
}

void ParallelParseBAM(const char *bam_input,
	vec<vec<pair<int, int> > > &lbps,
	vec<int> &scaffold_lens,
	vec<String> &scaffold_names,
	int thread_num = 8) {

	std::cerr << Date() << ": start to parse bam file" << std::endl;
	omp_set_num_threads(thread_num);
	samFile **bam_fp = new samFile*[thread_num];
	hts_idx_t **hts_idx = new hts_idx_t*[thread_num];
	bam1_t **bam_record = new bam1_t*[thread_num];

	// open bam file.
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < thread_num; ++i) {
		if ((bam_fp[i] = sam_open(bam_input, "r")) == nullptr) {
			std::cerr << "Can't open sam/bam file: " << bam_input << std::endl;
			exit(EXIT_FAILURE);
		}

		// get bam index content.
		if ((hts_idx[i] = sam_index_load(bam_fp[i], bam_input)) == nullptr) {
			std::cerr << "Can't get index for bam file" << std::endl;
			exit(EXIT_FAILURE);
		}
		bam_record[i] = bam_init1();
	}
	// get header info.
	bam_hdr_t *bam_header = nullptr;
	if ((bam_header = sam_hdr_read(bam_fp[0])) == nullptr) {
		std::cerr << Date() << ": Can't read bam header." << std::endl;
		Scram(1);
	}

	// clear lbp
	lbps.clear_and_resize(bam_header->n_targets);

	// record scaffold length and name
	scaffold_lens.clear();
	scaffold_names.clear();
	for (int i = 0; i < bam_header->n_targets; ++i) {
		scaffold_lens.push_back(bam_header->target_len[i]);
		scaffold_names.push_back(bam_header->target_name[i]);
	}

	int64_t read_cnt_global = 0, read_filtered_global = 0;
	int64_t low_qual_global = 0, dup_global = 0;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < bam_header->n_targets; ++i) {
		int thread_id = omp_get_thread_num();
		int64_t read_cnt = 0, read_filtered = 0, low_qual = 0, dup = 0;
		hts_itr_t *hts_itr = sam_itr_queryi(hts_idx[thread_id],
			i,
			0,
			bam_header->target_len[i]);
		while (sam_itr_next(bam_fp[thread_id],
			hts_itr,
			bam_record[thread_id]) >= 0) {
			++read_cnt;
			if (bam_record[thread_id]->core.qual < 60) {
				++low_qual;
				continue;
			}
			if (bam_record[thread_id]->core.flag &
				(0x400)) {
				++dup;
				continue;
			}
			if (bam_record[thread_id]->core.flag &
				(0x4 | 0x100 | 0x200 | 0x800)) {
				++read_filtered;
				continue;
			}
			String barcode_str =
				reinterpret_cast<char*>(bam_record[thread_id]->data);
			barcode_str = barcode_str.SafeAfterLast("#");
			barcode_t barcode_int = createBarcodeIndex(barcode_str);
			if (barcode_int == 0) {
				++read_filtered;
				continue;
			}
			// save barcode locations to lbp
			// should be thread-safe
			lbps[bam_record[thread_id]->core.tid].push(
				((bam_endpos(bam_record[thread_id]) +
					bam_record[thread_id]->core.pos) / 2),
				barcode_int);
		}
#pragma omp critical(read_count)
		{
			read_cnt_global += read_cnt;
			read_filtered_global += read_filtered;
			low_qual_global += low_qual;
			dup_global += dup;
		}
	}


	// sort lbp
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < lbps.isize(); ++i) {
		Sort(lbps[i]);
	}

	read_filtered_global += low_qual_global + dup_global;
	std::cerr << Date() << ": total reads: " << read_cnt_global
		<< std::endl;
	std::cerr << Date() << ": duplicated reads: " << dup_global
		<< std::endl;
	std::cerr << Date() << ": low quality: " << low_qual_global
		<< std::endl;
	std::cerr << Date() << ": filtered reads: " << read_filtered_global
		<< std::endl;
	std::cerr << Date() << ": qualified reads: "
		<< read_cnt_global - read_filtered_global << std::endl;
	std::cerr << Date() << ": end of parsing bam file" << std::endl;

	bam_hdr_destroy(bam_header);
	// free file handlers.
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < thread_num; ++i) {
		bam_destroy1(bam_record[i]);
		hts_idx_destroy(hts_idx[i]);
		sam_close(bam_fp[i]);
	}
}

void ParseFasta(const String &fasta_input,
	vec<String> &original_scaffolds,
	const vec<int> &scaffold_lens,
	const vec<String> &scaffold_names) {  // to check if length is correct

	std::cerr << Date() << ": start to parse fasta file" << std::endl;
	String command;
	if (fasta_input.EndsWith(".gz")) {
		command = "zcat " + fasta_input;
	}
	else {
		command = "cat " + fasta_input;
	}

	vec<String> scaffold_names_from_fasta;
	fast_pipe_ifstream input(command);
	String buf;
	buf.reserve(1024 * 1024);
	String scaffold;
	scaffold.reserve(*std::max_element(scaffold_lens.begin(),
		scaffold_lens.end()));
	size_t line = 0;
	while (1) {
		getline(input, buf);
		if (input.fail()) break;
		line++;
		// first or next seq

		if (buf.empty()) {
			continue;
		}
		if (buf.StartsWith(">")) {
			buf = buf.SafeAfter(">");
			buf = buf.SafeBefore(" ");
			scaffold_names_from_fasta.push_back(buf);
			if (scaffold.empty()) {
				continue;
			}
			// push current scaffold into vec
			original_scaffolds.push_back(scaffold);
			scaffold.clear();
		}
		else {
			scaffold += buf;
		}

	}
	// deal with last scaffold
	if (!scaffold.empty()) {
		original_scaffolds.push_back(scaffold);
		scaffold.clear();
	}

	// check if lengths are the same as in bam
	ForceAssertEq(original_scaffolds.size(), scaffold_lens.size());
	for (int i = 0; i < original_scaffolds.isize(); ++i) {
		if (original_scaffolds[i].isize() != scaffold_lens[i]) {
			std::cerr << Date() << ": #" << i << " scaffold's length not correct"
				<< std::endl;
		}
	}
	// check if names are the same as in bam
	ForceAssertEq(scaffold_names_from_fasta.size(), scaffold_names.size());
	for (int i = 0; i < scaffold_names_from_fasta.isize(); ++i) {
		if (scaffold_names_from_fasta[i] != scaffold_names[i]) {
			std::cerr << Date() << ": #" << i << " scaffold's name not correct"
				<< std::endl;
			std::cerr << Date() << scaffold_names_from_fasta[i] << ": "
				<< scaffold_names[i] << std::endl;
		}
	}
	std::cerr << Date() << ": scaffolds cnt = " << original_scaffolds.size()
		<< std::endl;
	std::cerr << Date() << ": end of parsing fasta file" << std::endl;
}

void GetGapArea(const vec<String> &original_scaffolds,
	vec<vec<pair<int, int> > > &gap_areas,
	int MIN_GAP) {

	std::cerr << Date() << ": start to create gap areas" << std::endl;
	int min_gap = 1000000000, max_gap = -1;
	gap_areas.clear_and_resize(original_scaffolds.size());
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < original_scaffolds.isize(); ++i) {
		const String &scaffold = original_scaffolds[i];
		vec<pair<int, int> > &gap_area = gap_areas[i];
		Bool inside_N_zone = False;
		int start, end;
		for (int j = 0; j < scaffold.isize(); ++j) {
			if (scaffold[j] == 'N' || scaffold[j] == 'n') {
				if (!inside_N_zone) {
					start = j;
					inside_N_zone = True;
				}
			}
			else {
				if (inside_N_zone) {
					end = j;
					ForceAssertLt(start, end);
					int gap_len = end - start;
					if (gap_len < min_gap) {
						min_gap = gap_len;
					}
					if (gap_len > max_gap) {
						max_gap = gap_len;
					}

					if (gap_len >= MIN_GAP) {
						gap_area.push(start, end);
					}
					inside_N_zone = False;
				}
			}
		}
		// if a scaffold ends with Ns
		if (inside_N_zone) {
			end = scaffold.isize();
			ForceAssertLt(start, end);
			int gap_len = end - start;
			if (gap_len < min_gap) {
				min_gap = gap_len;
			}
			if (gap_len > max_gap) {
				max_gap = gap_len;
			}

			if (gap_len >= MIN_GAP) {
				gap_area.push(start, end);
			}
		}
	}

	std::cerr << Date() << ": min_gap = " << min_gap << std::endl;
	std::cerr << Date() << ": max_gap = " << max_gap << std::endl;
	std::cerr << Date() << ": end of creating gap areas" << std::endl;
}


void BucketLines(const vec<int>& llens,
	vec<vec<int> >& buckets,
	const int min_len)
{
	buckets.clear();
	const int bsize = 1000000;
	int nd = llens.size();
	vec<int> llensx(llens), ids(nd, vec<int>::IDENTITY);
	ParallelReverseSortSync(llensx, ids);
	for (int i = 0; i < nd; i++)
	{
		if (llensx[i] < min_len) continue;
		int64_t sum = llensx[i];
		int j;
		for (j = i + 1; j < nd; j++)
		{
			if (llensx[j] < min_len) continue;
			if (sum >= bsize) break;
			sum += llensx[j];
		}
		vec<int> b;
		for (int k = i; k < j; k++)
		{
			if (llensx[k] < min_len) continue;
			b.push_back(ids[k]);
		}
		buckets.push_back(b);
		i = j - 1;
	}
}

void LbpOnReverseComplement(vec<pair<int, int> > &lbp, int len) {
	for (int i = 0; i < lbp.isize(); ++i) {
		lbp[i].first = len - lbp[i].first - 1;
	}
	Sort(lbp);
}

void GapAreaOnReverseComplement(vec<pair<int, int> > &gap_area, int len) {
	for (int i = 0; i < gap_area.isize(); ++i) {
		gap_area[i].first = len - gap_area[i].first - 1;
		gap_area[i].second = len - gap_area[i].second - 1;
	}
	Sort(gap_area);
}

// create rc gap_area and lbpx for linked scaffolds
void GapEstPrepare(const vec<vec<int> > &linked_scaffolds,
	const vec<vec<pair<int, int> > > &gap_areas,
	const vec<vec<pair<int, int> > > &lbps,
	const vec<int> &scaffold_lens,
	vec<vec<int> > &bc_gap_pos,
	vec<vec<pair<int, int> > > &gap_areas_dummy,
	vec<vec<pair<int, int> > > &lbpxs) {

	lbpxs.clear_and_resize(linked_scaffolds.size());
	bc_gap_pos.clear_and_resize(linked_scaffolds.size());
	gap_areas_dummy.clear_and_resize(linked_scaffolds.size());
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < linked_scaffolds.isize(); ++i) {
		const vec<int> &linked_scaffold = linked_scaffolds[i];
		if (linked_scaffold.solo()) {
			continue;
		}
		int len_cur = 0;
		for (int j = 0; j < linked_scaffold.isize(); ++j) {
			int sid = linked_scaffold[j];
			vec<pair<int, int> > lbp_tmp = lbps[sid / 2];
			vec<pair<int, int> > gap_area_tmp = gap_areas[sid / 2];
			if (sid % 2 != 0) { // reverese complement
				LbpOnReverseComplement(lbp_tmp, scaffold_lens[sid / 2]);
				GapAreaOnReverseComplement(gap_area_tmp,
					scaffold_lens[sid / 2]);
			}
			for (int k = 0; k < lbp_tmp.isize(); ++k) {
				lbp_tmp[k].first += len_cur;
			}
			for (int k = 0; k < gap_area_tmp.isize(); ++k) {
				gap_area_tmp[k].first += len_cur;
				gap_area_tmp[k].second += len_cur;
			}
			lbpxs[i].append(lbp_tmp);
			gap_areas_dummy[i].append(gap_area_tmp);
			len_cur += scaffold_lens[sid / 2];
			bc_gap_pos[i].push_back(len_cur); // last ele is length of l_scaff

		}
		// just in case, actually no need to sort
		Sort(lbpxs[i]);
		Sort(gap_areas_dummy[i]);
	}
	/*
	for (int i = 0; i < gap_areas_dummy.isize(); ++i) {
		if (gap_areas_dummy[i].empty()) continue;
		for (int j = 0; j < gap_areas_dummy[i].isize(); ++j) {
			std::cerr << gap_areas_dummy[i][j].first << ':' << gap_areas_dummy[i][j].second << ' ';
		}
		std::cerr << std::endl;
	}
	*/

}

void GapEstTest(const vec<vec<pair<int, int> > > &gap_areas,
	const vec<int> &scaffold_lens,
	const vec<vec<pair<int, int> > > &lbps,
	const vec<pair<int, double> > &gb) {

	std::cerr << Date() << ": start to test gap distance est" << std::endl;
	const int WINDOW = 10000;
	const int MULT = 50;
	const int GAP_DELTA = 5000;
	const int MIN_GAP = 1000;
	const int MIN_POINTS = 2;

	vec<vec<int> > bc_gap_distances(scaffold_lens.size(), vec<int>());
	for (int i = 0; i < bc_gap_distances.isize(); ++i) {
		vec<int> &bc_gap_distance = bc_gap_distances[i];
		bc_gap_distance.assign(gap_areas[i].size(), -1);
	}

#pragma omp parallel for schedule(dynamic,1)
	for (int l = 0; l < gap_areas.isize(); l++) {
		const vec<pair<int, int> > &gap_area = gap_areas[l];
		if (gap_area.empty()) {
			continue;
		}

		const vec< pair<int, int> > &pb = lbps[l];
		vec<int> lr;
		// last element of bc_gap_pos[id] is the length of linked scaff
		for (int pos_id = 0; pos_id < gap_area.isize(); ++pos_id) {
			pair<int, int> pos = gap_area[pos_id];
			int start = pos.first, end = pos.second;
			if (end - start < MIN_GAP) {
				continue;
			}
			vec<int> bcount, lefts, rights;
			vec<double> fracs;
			int left1 = start - WINDOW, left2 = start;
			int right1 = end, right2 = end + WINDOW;
			if (left1 < 0 || right2 > scaffold_lens[l]) continue; // FAIL
			Bool at_gap = False;
			for (auto &area : gap_areas[l]) {
				if (area.first == start && area.second == end) {
					continue;
				}
				int l_bound = max(area.first, left1);
				int r_bound = min(area.second, right2);
				if (l_bound - r_bound < 0) { // overlapped
					at_gap = True;
					break;
				}
			}
			if (at_gap) continue; // FAIL
			int low = LowerBound1(pb, left1);
			for (int j = low; j < pb.isize(); j++)
			{
				int p = pb[j].first;
				if (p > right2) break;
				int b = pb[j].second;
				if (p >= left1 && p < left2) lefts.push_back(b);
				if (p >= right1 && p < right2) rights.push_back(b);
			}

			Sort(lefts), Sort(rights);
			vec<Bool> del_left(lefts.size(), False);
			for (int j = 0; j < lefts.isize(); j++)
			{
				int k;
				for (k = j + 1; k < lefts.isize(); k++)
					if (lefts[k] != lefts[j]) break;
				if (k - j < MIN_POINTS)
					for (int l = j; l < k; l++) del_left[l] = True;
				j = k - 1;
			}
			EraseIf(lefts, del_left);
			vec<Bool> del_right(rights.size(), False);
			for (int j = 0; j < rights.isize(); j++)
			{
				int k;
				for (k = j + 1; k < rights.isize(); k++)
					if (rights[k] != rights[j]) break;
				if (k - j < MIN_POINTS)
					for (int l = j; l < k; l++) del_right[l] = True;
				j = k - 1;
			}
			EraseIf(rights, del_right);

			UniqueSort(lefts), UniqueSort(rights);
			lr = lefts;
			lr.append(rights);
			UniqueSort(lr);
			int bridges = MeetSize(lefts, rights);
			if (lr.empty()) {
				continue;
			}
			double bridge_frac = double(bridges) / lr.size();

			// If linking is too weak, don't set gap.  In some cases these
			// instances would be due to misassemblies.

			if (gb.nonempty() && bridge_frac < gb.back().second / 2) {
				continue;
			}

			int best_gap = -1;
			double score = 1000000000;
			for (int i = 0; i < gb.isize(); i++)
			{
				int gap = gb[i].first;
				double count = gb[i].second;
				if (Abs(count - bridge_frac) < score)
				{
					score = Abs(count - bridge_frac);
					best_gap = gap;
				}
			}
			if (best_gap < 0) continue;
			best_gap = Max(best_gap, MIN_GAP);
			// Set the gap
			bc_gap_distances[l][pos_id] = best_gap;
		}
	}

	std::cerr << "gap_actual\tgap_est" << std::endl;
	int perfect = 0, tiny_diff = 0, too_much = 0;

	for (int i = 0; i < bc_gap_distances.isize(); ++i) {
		for (int j = 0; j < bc_gap_distances[i].isize(); ++j) {
			int start = gap_areas[i][j].first;
			int end = gap_areas[i][j].second;
			int actual_gap = end - start;
			if (actual_gap < MIN_GAP) continue;
			std::cerr << actual_gap << '\t'
				<< bc_gap_distances[i][j] << std::endl;
			if (bc_gap_distances[i][j] == -1 && actual_gap == 3000) {
				++perfect;
				continue;
			}
			int gap_diff = Abs(actual_gap - bc_gap_distances[i][j]);
			if (gap_diff == 0) {
				++perfect;
			}
			else if (gap_diff <= 5000) {
				++tiny_diff;
			}
			else {
				++too_much;
			}
		}
	}
	double total = perfect + tiny_diff + too_much;
	std::cerr << "total: " << total << std::endl;
	std::cerr << "perfect: " << perfect / total << std::endl;
	std::cerr << "tiny: " << tiny_diff / total << std::endl;
	std::cerr << "too much: " << too_much / total << std::endl;
	std::cerr << Date() << ": end of testing gap distance est" << std::endl;

}

void GaprikaExternal(vec<vec<int> > &linked_scaffolds,
	const vec<vec<pair<int, int> > > &gap_areas,
	const vec<int> &scaffold_lens,  // original scaffold length, without gap
	const vec<vec<pair<int, int> > > &lbps,
	vec<vec<int> > &bc_gap_distances,
	const int round,
	const int MAX_GAP,
	const int VERBOSITY, int window)
{
	std::cerr << Date() << ": start to gaprika" << std::endl;
	// Heuristics.

	int WINDOW = window;
	const int MULT = 50;
	const int GAP_DELTA = 5000;
	const int MIN_GAP = 400;
	const int MIN_POINTS = 2;

	bc_gap_distances.assign(linked_scaffolds.size(), vec<int>());
	for (int i = 0; i < bc_gap_distances.isize(); ++i) {
		vec<int> &bc_gap_distance = bc_gap_distances[i];
		bc_gap_distance.assign(linked_scaffolds[i].size() - 1, -1);
	}
	double clock = WallClockTime();

	// Range over gap sizes.

	vec<vec<int> > buckets;
	BucketLines(scaffold_lens, buckets, 2 * WINDOW);
	cerr << Date() << ": starting loop" << endl;
	if (VERBOSITY >= 1) cerr << endl;
	vec<pair<int, double> > gb;
	for (int gap = 0; gap <= MAX_GAP; gap += GAP_DELTA)
	{
		vec<int> nbridges;
		vec<double> fracs;
#pragma omp parallel for schedule(dynamic,1)
		for (int bu = 0; bu < buckets.isize(); bu++)
		{
			vec<int> bcount, lefts, rights, lr;
			vec<double> fr;
			for (auto l : buckets[bu])
			{
				if (scaffold_lens[l] < 2 * WINDOW + gap) continue;
				const vec< pair<int, int> > &pb = lbps[l];
				bcount.clear();
				fr.clear();
				// Sort(pb);
				for (int i = WINDOW; i <= scaffold_lens[l] - WINDOW - gap;
					i += WINDOW * MULT)
				{
					int left1 = i - WINDOW, left2 = i;
					int right1 = i + gap, right2 = i + WINDOW + gap;
					Bool at_gap = False;
					for (auto &area : gap_areas[l]) {
						int l_bound = max(area.first, left1);
						int r_bound = min(area.second, right2);
						// if ( b >= left1 && b <= right2 ) at_gap = True;
						if (l_bound - r_bound < 0) { // overlapped
							at_gap = True;
							break;
						}
					}
					if (at_gap) continue;
					int low = LowerBound1(pb, left1);
					lefts.clear(), rights.clear();
					for (int j = low; j < pb.isize(); j++)
					{
						int p = pb[j].first;
						if (p > right2) break;
						int b = pb[j].second;
						if (p >= left1 && p < left2) lefts.push_back(b);
						if (p >= right1 && p < right2)
							rights.push_back(b);
					}

					Sort(lefts), Sort(rights);
					vec<Bool> del_left(lefts.size(), False);
					for (int j = 0; j < lefts.isize(); j++)
					{
						int k;
						for (k = j + 1; k < lefts.isize(); k++)
							if (lefts[k] != lefts[j]) break;
						if (k - j < MIN_POINTS)
							for (int l = j; l < k; l++) del_left[l] = True;
						j = k - 1;
					}
					EraseIf(lefts, del_left);
					vec<Bool> del_right(rights.size(), False);
					for (int j = 0; j < rights.isize(); j++)
					{
						int k;
						for (k = j + 1; k < rights.isize(); k++)
							if (rights[k] != rights[j]) break;
						if (k - j < MIN_POINTS)
							for (int l = j; l < k; l++) del_right[l] = True;
						j = k - 1;
					}
					EraseIf(rights, del_right);

					UniqueSort(lefts), UniqueSort(rights);
					lr = lefts;
					lr.append(rights);
					UniqueSort(lr);
					int bridges = MeetSize(lefts, rights);
					if (lr.empty()) continue;
					double bridge_frac = double(bridges) / lr.size();
					bcount.push_back(bridges);
					fr.push_back(bridge_frac);
				}
#pragma omp critical(append_bridge)
				{
					nbridges.append(bcount);
					fracs.append(fr);
				}
			}
		}
		Sort(fracs);
		int count = fracs.size();
		if (count > 0)
		{
			double m = Mean(fracs);
			double dev = StdDev(fracs, m);
			if (VERBOSITY >= 1) PRINT5(gap, count, m, dev, Median(fracs));
			gb.push(gap, m);
		}
	}
	if (VERBOSITY >= 1) cerr << endl;
	cerr << Date() << ": done, time used = " << TimeSince(clock) << endl;

	// Now estimate sizes of bc-only gaps.
	vec<vec<int> > bc_gap_pos;
	vec<vec<pair<int, int> > > gap_areas_dummy;
	vec<vec<pair<int, int> > > lbpxs;
	GapEstPrepare(linked_scaffolds, gap_areas, lbps, scaffold_lens,
		bc_gap_pos, gap_areas_dummy, lbpxs);

	cerr << Date() << ": estimating gap sizes" << endl;
	if (VERBOSITY >= 1) cerr << endl;
	// #pragma omp parallel for schedule(dynamic,1)
	for (int l = 0; l < linked_scaffolds.isize(); l++) {

		const vec<int> &linked_scaffold = linked_scaffolds[l];
		if (linked_scaffold.solo()) {
			continue;
		}
		/*
		vec< pair<int,int>> pb;
		for ( auto x : lbpx[l] ) pb.push( x.second, x.first );
		Sort(pb);
		*/
		const vec< pair<int, int> > &pb = lbpxs[l];
		vec<int> lr;
		// last element of bc_gap_pos[id] is the length of linked scaff
		for (int pos_id = 0; pos_id < bc_gap_pos[l].isize() - 1; ++pos_id) {
			int pos = bc_gap_pos[l][pos_id];
			vec<int> gap_len_est(round, -1);
			vec<Bool> at_gap_global(round, False);
			for (int pass = 1; pass <= round; ++pass) {

				vec<int> bcount, lefts, rights;
				vec<double> fracs;
				int linked_scaffold_len = bc_gap_pos[l].back();
				int left1 = pos - pass * WINDOW, left2 = pos - (pass - 1) * WINDOW;
				int right1 = pos + (pass - 1) * WINDOW, right2 = pos + pass * WINDOW;
				if (left1 < 0 || right2 > linked_scaffold_len) continue; // FAIL
				Bool at_gap = False;
				for (auto &area : gap_areas_dummy[l]) {
					int l_bound = max(area.first, left1);
					int r_bound = min(area.second, right2);
					if (l_bound - r_bound < 0) { // overlapped
						at_gap = True;
						break;
					}
				}
				if (at_gap) {
					at_gap_global[pass - 1] = True;
					break; // FAIL
				}
				int low = LowerBound1(pb, left1);
				for (int j = low; j < pb.isize(); j++)
				{
					int p = pb[j].first;
					if (p > right2) break;
					int b = pb[j].second;
					if (p >= left1 && p < left2) lefts.push_back(b);
					if (p >= right1 && p < right2) rights.push_back(b);
				}

				Sort(lefts), Sort(rights);
				vec<Bool> del_left(lefts.size(), False);
				for (int j = 0; j < lefts.isize(); j++)
				{
					int k;
					for (k = j + 1; k < lefts.isize(); k++)
						if (lefts[k] != lefts[j]) break;
					if (k - j < MIN_POINTS)
						for (int l = j; l < k; l++) del_left[l] = True;
					j = k - 1;
				}
				EraseIf(lefts, del_left);
				vec<Bool> del_right(rights.size(), False);
				for (int j = 0; j < rights.isize(); j++)
				{
					int k;
					for (k = j + 1; k < rights.isize(); k++)
						if (rights[k] != rights[j]) break;
					if (k - j < MIN_POINTS)
						for (int l = j; l < k; l++) del_right[l] = True;
					j = k - 1;
				}
				EraseIf(rights, del_right);

				UniqueSort(lefts), UniqueSort(rights);
				lr = lefts;
				lr.append(rights);
				UniqueSort(lr);
				int bridges = MeetSize(lefts, rights);
				if (lr.empty()) {
					continue;
				}
				double bridge_frac = double(bridges) / lr.size();

				// If linking is too weak, don't set gap.  In some cases these
				// instances would be due to misassemblies.

				if (gb.nonempty() && bridge_frac < gb.back().second / 2) {
					continue;
				}

				int best_gap = -1;
				double score = 1000000000;
				for (int i = 0; i < gb.isize(); i++)
				{
					int gap = gb[i].first;
					double count = gb[i].second;
					if (Abs(count - bridge_frac) < score)
					{
						score = Abs(count - bridge_frac);
						best_gap = gap;
					}
				}
				if (best_gap < 0) continue;
				best_gap = Max(best_gap, MIN_GAP);
				// Set the gap candidate
				gap_len_est[pass - 1] = best_gap;
				// bc_gap_distances[l][pos_id] = best_gap;
			}

			if (VERBOSITY >= 1) {
				for (int pass = 1; pass <= round; ++pass) {
					//std::cerr << int(at_gap_global[pass - 1]) << ' ';
				}
				std::cerr << std::endl;
				for (int pass = 1; pass <= round; ++pass) {
					int gap_cur = gap_len_est[pass - 1] - 2 * (pass - 1) * WINDOW;
					//std::cerr << gap_cur << ' ';
				}
				//std::cerr << std::endl;
			}

			int best_gap = 0, times = 0;
			for (int pass = 1; pass <= round; ++pass) {
				if (at_gap_global[pass - 1]) {
					break;
				}
				if (gap_len_est[pass - 1] == -1) {
					continue;
				}
				int gap_cur = gap_len_est[pass - 1] - 2 * (pass - 1) * WINDOW;
				if (gap_cur >= 0) {
					best_gap += gap_cur;
					++times;
				}
			}
			if (times == 0) continue;
			bc_gap_distances[l][pos_id] = best_gap / times;
		}
	}

	// deal with "-1" gaps, split seq between it
#if 1 
	for (int i = 0; i < bc_gap_distances.isize(); ++i) {
		int prev = 0;
		for (int j = 0; j < bc_gap_distances[i].isize(); ++j) {
			if (bc_gap_distances[i][j] == -1) {
				bc_gap_distances[i][j] = 30000;
			}
		}
	}
#else 

	vec<vec<int> > linked_scaffolds_tmp;
	vec<vec<int> > bc_gap_distances_tmp;
	for (int i = 0; i < bc_gap_distances.isize(); ++i) {
		int prev = 0;
		for (int j = 0; j < bc_gap_distances[i].isize(); ++j) {
			if (bc_gap_distances[i][j] == -1) {
				linked_scaffolds_tmp.push_back(
					vec<int>(linked_scaffolds[i].begin() + prev,
						linked_scaffolds[i].begin() + j + 1));
				bc_gap_distances_tmp.push_back(
					vec<int>(bc_gap_distances[i].begin() + prev,
						bc_gap_distances[i].begin() + j));
				prev = j + 1;

			}
		}
		// deal with last part
		linked_scaffolds_tmp.push_back(
			vec<int>(linked_scaffolds[i].begin() + prev,
				linked_scaffolds[i].end()));
		bc_gap_distances_tmp.push_back(
			vec<int>(bc_gap_distances[i].begin() + prev,
				bc_gap_distances[i].end()));

	}
	bc_gap_distances.swap(bc_gap_distances_tmp);
	linked_scaffolds.swap(linked_scaffolds_tmp);
	// GapEstTest(gap_areas, scaffold_lens, lbps, gb);
#endif

	std::cerr << Date() << ": end of gaprika" << std::endl;
}

void OutoutSeq(std::ofstream &f_out,
	String &seq) {
	int cnt = 0;
	for (int i = 0; i < seq.isize(); ++i) {
		f_out << seq[i];
		++cnt;
		if (cnt % 80 == 0) {
			f_out << std::endl;
		}
	}
	if (seq.size() % 80 != 0) {
		f_out << std::endl;
	}
}

char RCBase(char base) {
	if (base == 'A' || base == 'a') {
		return 'T';
	}
	else if (base == 'C' || base == 'c') {
		return 'G';
	}
	else if (base == 'T' || base == 't') {
		return 'A';
	}
	else if (base == 'G' || base == 'g') {
		return 'C';
	}
	return base;
}

String ReverseComplementSeq(const String &ori) {
	String rc;
	rc.reserve(ori.size());
	for (auto it = ori.rbegin(); it != ori.rend(); ++it) {
		// std::cerr << "single base: " << *it << std::endl;
		// std::cerr << "reverse base: " << RCBase(*it) << std::endl;
		rc.push_back(RCBase(*it));
	}
	return rc;
}

void TXT2GZ(const String& fn) {
	ifstream in;
	std::string data_fn;
	in.open(fn.c_str(), std::ios::in);
	if (in) {
		if (!in.seekg(0, std::ios::end)) {
			cerr << "Error seeking in file " << fn << endl;
			Scram(1);
		}
		data_fn.reserve(in.tellg());
		if (!in.seekg(0, std::ios::beg)) {
			cerr << "Error seeking in file " << fn << endl;
			Scram(1);
		}
		data_fn.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
		in.close();
	}
	else {
		cerr << "Error opening " << fn << " for reading" << endl;
		Scram(1);
	}
	FstreamWriteGZ(fn, data_fn.c_str());
	std::cerr << Date() << ": finish writing gzipped object" << endl;
}

void OutputLinkedScaffold(const String &fasta_out,
	const vec<vec<int> > &linked_scaffolds,
	const vec<String> &scaffold_names,
	vec<String> &original_scaffolds,
	const vec<vec<int> > &bc_gap_distances,
	Bool ZIP) {
	std::ofstream f_out;
	f_out.open(fasta_out.c_str(), std::ofstream::out);
	if (!f_out.is_open()) {
		std::cerr << Date() << ": can't open output fasta file" << std::endl;
		Scram(1);
	}
	String seq;
	seq.reserve(1024 * 1024 * 1024);
	for (int i = 0; i < linked_scaffolds.isize(); ++i) {
		const vec<int> &linked_scaffold = linked_scaffolds[i];
		const vec<int> & bc_gap_distance = bc_gap_distances[i];
		f_out << '>' << i << ' ';
		if (linked_scaffold.solo()) {
			int id = linked_scaffold.front();
			f_out << scaffold_names[id / 2] << std::endl;
			// no rc case
			OutoutSeq(f_out, original_scaffolds[id / 2]);
			continue;
		}
		int id = linked_scaffold.front();
		f_out << scaffold_names[id / 2];
		if (id % 2 != 0) {
			seq.append(ReverseComplementSeq(original_scaffolds[id / 2]));
		}
		else {
			seq.append(original_scaffolds[id / 2]);
		}
		for (int j = 1; j < linked_scaffold.isize(); ++j) {
			int gap_cnt = bc_gap_distance[j - 1];
			ForceAssertNe(gap_cnt, -1);
			seq.append(gap_cnt, 'N');
			id = linked_scaffold[j];
			f_out << " + " << gap_cnt << " N + " << scaffold_names[id / 2];
			if (id % 2 != 0) {
				seq.append(ReverseComplementSeq(original_scaffolds[id / 2]));
			}
			else {
				seq.append(original_scaffolds[id / 2]);
			}
		}
		f_out << std::endl;
		OutoutSeq(f_out, seq);
		seq.clear();
	}
	f_out.close();
	if (ZIP) TXT2GZ(fasta_out);
}

int main(int argc, char **argv) {

	RunTime();
	double clock = WallClockTime();

	BeginCommandArguments;
	CommandArgument_Bool_OrDefault_Doc(TRACK_SOME_MEMORY, False,
		"track some memory allocations, for debugging");
	CommandArgument_String_OrDefault_Doc(DIR, ".", "input directory");

	CommandArgument_String_OrDefault_Doc(BAM_IN, "", "bam input");
	CommandArgument_String_OrDefault_Doc(FASTA_IN, "", "fasta input");
	CommandArgument_String_OrDefault_Doc(PAF_IN, "", "minimap2 paf input");
	CommandArgument_String_OrDefault_Doc(LWML_IN, "", "lwml input");
	CommandArgument_String_OrDefault_Doc(FASTA_OUT, "", "fasta output");
	CommandArgument_Int_OrDefault_Doc(EST_ROUND, 3, "how many rounds to est gap");
	CommandArgument_Int_OrDefault_Doc(MIN_GAP, 400,
		"min gap for gather gap info");
	CommandArgument_Int_OrDefault_Doc(MAX_GAP, 100000, "max gap for GAPRIKA");
	CommandArgument_Bool_OrDefault_Doc(ZIP, True, "output in gzip format");

	CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
		"Number of threads.  By default, the number of processors online.");
	CommandArgument_Double_OrDefault_Doc(MAX_MEM_GB, 0,
		"if specified, maximum *suggested* RAM use in GB; in some cases may be "
		"exceeded by our code");

	EndCommandArguments;

	if (TRACK_SOME_MEMORY) DeclareThatWeAreTrackingSomeMemory();
	// Set computational limits, etc.

	SetThreads(NUM_THREADS, False);
	SetMaxMemoryGBCheck(MAX_MEM_GB);

	/*
	// random test
	String seq = "ACTGNNNTGAC";
	std::cerr << ReverseComplementSeq(seq) << std::endl;
	Scram(0);
	*/

	// parse all kinds of input
	vec<vec<pair<int, int> > > lbps;
	vec<int> scaffold_lens;
	vec<String> scaffold_names;
	int window=20000;
	if (LWML_IN != "")
	{
		BinaryReader::readFile(LWML_IN.c_str(), &window);
	}
	

	window /= 4;
	//cout << window << endl;
	// ParseBAM(BAM_IN.c_str(), lbps, scaffold_lens, scaffold_names);
	ParallelParseBAM(BAM_IN.c_str(), lbps, scaffold_lens, scaffold_names,
		NUM_THREADS);

	vec<String> original_scaffolds;
	ParseFasta(FASTA_IN, original_scaffolds, scaffold_lens, scaffold_names);

	vec<vec<pair<int, int> > > gap_areas;
	GetGapArea(original_scaffolds, gap_areas, MIN_GAP);

	vec<Bool> non_dup;
	dr_wu::parse_paf(PAF_IN.c_str(), scaffold_names, non_dup);

	vec<vec<int> > linked_scaffolds;
	dr_wu::linear_scaffolding(lbps, scaffold_lens, gap_areas, non_dup, scaffold_names, linked_scaffolds, window);

	vec<vec<int> > bc_gap_distances;
	int VERBOSITY = 1;
	GaprikaExternal(linked_scaffolds, gap_areas, scaffold_lens, lbps,
		bc_gap_distances, EST_ROUND, MAX_GAP, VERBOSITY, window);


	OutputLinkedScaffold(FASTA_OUT,
		linked_scaffolds,
		scaffold_names,
		original_scaffolds,
		bc_gap_distances,
		ZIP);

	Scram(0);
}

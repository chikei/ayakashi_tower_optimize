#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <fstream>
#include "comb.hpp"
#include <cstdint>
#include <iostream>
#include <chrono>
#include <boost/tokenizer.hpp>
#include <utility>

const std::string phantom = u8"\u5996";
const std::string divina = u8"\u795E";
const std::string anima = u8"\u4E5D\u5341\u4E5D";
const std::string all = u8"\u5168";
const std::string self_buff = u8"\u81EA\u8EAB";
const std::string team_buff = u8"\u6211\u65B9";
const std::string enemy_buff = u8"\u6575\u65B9";
const std::string buff_atk_def = u8"\u653B\u9632";
const std::string buff_atk = u8"\u653B";
const std::string buff_def = u8"\u9632";
const std::vector<std::vector<bool>> selections{ { 0, 0, 0, 0, 0}, 
  { 1, 0, 0, 0, 0},
  { 0, 1, 0, 0, 0},
  { 1, 1, 0, 0, 0},
  { 0, 0, 1, 0, 0},
  { 1, 0, 1, 0, 0},
  { 0, 1, 1, 0, 0},
  { 1, 1, 1, 0, 0},
  { 0, 0, 0, 1, 0},
  { 1, 0, 0, 1, 0},
  { 0, 1, 0, 1, 0},
  { 1, 1, 0, 1, 0},
  { 0, 0, 1, 1, 0},
  { 1, 0, 1, 1, 0},
  { 0, 1, 1, 1, 0},
  { 1, 1, 1, 1, 0},
  { 0, 0, 0, 0, 1},
  { 1, 0, 0, 0, 1},
  { 0, 1, 0, 0, 1},
  { 1, 1, 0, 0, 1},
  { 0, 0, 1, 0, 1},
  { 1, 0, 1, 0, 1},
  { 0, 1, 1, 0, 1},
  { 1, 1, 1, 0, 1},
  { 0, 0, 0, 1, 1},
  { 1, 0, 0, 1, 1},
  { 0, 1, 0, 1, 1},
  { 1, 1, 0, 1, 1},
  { 0, 0, 1, 1, 1},
  { 1, 0, 1, 1, 1},
  { 0, 1, 1, 1, 1},
  { 1, 1, 1, 1, 1} };

struct monster{
  int id;
  std::string name;
  std::string lv;
  int cost;
  std::string attr;
  int attack;
  int skill_lv = 0;
  int skill_base = 0;
  std::string skill_target = "";
  std::string skill_attr = "";
  std::string skill_buff = "";
  int skill_freq = 50;
};

struct enemy{
  std::string attr;
  std::string name;
  int defense;
  int extra_defense = 0;
};

struct enemy_reduce_attack{
  std::string attr;
  std::string name;
  int defense;
  int skill;
  std::string skill_attr;
};

std::vector<std::string> split(const std::string& line){
  boost::char_separator<char> sep(",", "", boost::keep_empty_tokens);
  boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
  std::vector<std::string> data(tokens.begin(), tokens.end());
  return data;
}

std::shared_ptr<monster> parse_monster(const std::string& line){
  std::shared_ptr<monster> nptr;
  auto data = split(line);
  auto m = std::make_shared<monster>();
  if(data.size() < 11) return nptr;
  m->id = std::stoi(data[1]);
  m->name = data[2];
  m->lv = data[3];
  m->cost = std::stoi(data[4]);
  m->attr = data[5];
  m->attack = std::stoi(data[7]);

  if(!data[9].empty())
    m->skill_lv = std::stoi(data[9]);

  if(!data[10].empty()){
    if(data[10].find(u8"\u5C0F\u5E45") != std::string::npos){
      m->skill_base = 2;
    } else if(data[10].find(u8"\u5927\u5E45") != std::string::npos){
      m->skill_base = 9;
    } else if(data[10].find(u8"\u7279\u5927") != std::string::npos){
      m->skill_base = 19;
    } else {
      m->skill_base = 7;
    }
  }

  if(data[10].find(u8"\u81EA\u8EAB") != std::string::npos){
    m->skill_target = self_buff;
  } else if(data[10].find(u8"\u6575\u65B9") != std::string::npos){
    m->skill_target = enemy_buff;
  } else if(data[10].find(u8"\u63D0\u9AD8") != std::string::npos){
    m->skill_target = team_buff;
  }

  if(data[10].find(u8"\u5996\u9B54") != std::string::npos){
    m->skill_attr = phantom;
  } else if(data[10].find(u8"\u795E\u9748") != std::string::npos){
    m->skill_attr = divina;
  } else if(data[10].find(u8"\u4E5D\u5341\u4E5D") != std::string::npos){
    m->skill_attr = anima;
  }

  if(data[10].find(u8"\u653B\u64CA\u529B") != std::string::npos && 
      data[10].find(u8"\u9632\u79A6\u529B") != std::string::npos){
    m->skill_buff = buff_atk_def;
  } else if(data[10].find(u8"\u653B\u64CA\u529B") != std::string::npos){
    m->skill_buff = buff_atk;
  } else if(data[10].find(u8"\u9632\u79A6\u529B") != std::string::npos){
    m->skill_buff = buff_def;
  }

  if(data[10].find(u8"\u9AD8\u767C\u52D5\u7387") != std::string::npos){
    m->skill_freq = 70;
  }

  return m;
}

std::shared_ptr<enemy> parse_enemy(const std::string& line){
  std::shared_ptr<enemy> nptr;
  auto data = split(line);
  auto e = std::make_shared<enemy>();
  if(data.size() < 3)
    return nptr;
  e->attr = data[0];
  e->defense = std::stoi(data[2]);
  e->name = data[1];
  if(data.size() > 3)
    e->extra_defense = std::stoi(data[3]);
  return e;
}

std::shared_ptr<enemy_reduce_attack> parse_enemy_reduce_attack(const std::string& line){
  std::shared_ptr<enemy_reduce_attack> nptr;
  auto data = split(line);
  auto e = std::make_shared<enemy_reduce_attack>();
  if(data.size() < 5)
    return nptr;
  e->attr = data[0];
  e->defense = std::stoi(data[2]);
  e->name = data[1];
  e->skill = std::stoi(data[3]);
  e->skill_attr = data[4];
  return e;
}

double attr_advantage_bonus(int base, bool is_adv){
  return is_adv ? base * 0.15 : 0.0;
}

bool is_attr_advantage(const std::shared_ptr<monster> m, const std::string& attr){
  return (m->attr == phantom && attr == anima) ||
    (m->attr == anima && attr == divina) ||
    (m->attr == divina && attr == phantom);
}

double attr_same_bonus(const std::vector<std::shared_ptr<monster>>& our,
    std::shared_ptr<monster> m, bool is_adv){
  int count = 0;
  for(const auto& c: our) if(c->attr == m->attr) ++count;

  if(count <= 3) return 0.0;
  else if(count == 4){
    double rate = is_adv ? 0.02867 : 0.024999;
    double bonus = 0.0;
    for(const auto& o: our){
      if(o->attr == m->attr)
        bonus += o->attack * rate;
    }
    return bonus;
  } else {
    double rate = is_adv ? 0.05817 : 0.05056;
    double bonus = 0.0;
    for(const auto& o: our) bonus += o->attack * rate;
    return bonus;
  }
}

double compute_basic_atk(const std::vector<std::shared_ptr<monster>>& our,
    const std::shared_ptr<enemy> target,
    double base, int leader){
  auto is_adv = is_attr_advantage(our[leader], target->attr);
  return base + attr_advantage_bonus(base, is_adv) + 
    attr_same_bonus(our, our[leader], is_adv) - target->extra_defense;
}

template <typename E>
int try_attack(const std::vector<std::shared_ptr<monster>>& our,
    const std::shared_ptr<E> target){
  int base = 0;
  for(const auto& m : our) { base += m->attack; }
  for(int i = 0; i < 5; i++){
    auto atk = compute_basic_atk(our, target, base, i);
    if(atk > target->defense) return i;
  }
  return -1;
}

double compute_reduce(const std::vector<std::shared_ptr<monster>>& our, int reduce, 
    const std::string& reduced_attr){
  double r = 0.0;
  for(const auto& m: our){
    if(m->attr == reduced_attr || reduced_attr == all)
      r += m->attack * reduce / 100.0;
  }
  return r;
}

double compute_basic_atk(const std::vector<std::shared_ptr<monster>>& our,
    const std::shared_ptr<enemy_reduce_attack> target,
    double base, int leader){
  auto is_adv = is_attr_advantage(our[leader], target->attr);
  return base + attr_advantage_bonus(base, is_adv) + 
    attr_same_bonus(our, our[leader], is_adv) -
    compute_reduce(our, target->skill, target->skill_attr);
}

double compute_skill(const std::vector<std::shared_ptr<monster>>& our, 
    std::shared_ptr<monster> m, int defense, const std::string& attr){
  if(m->skill_lv == 0) return 0.0;
  else {
    auto str = (m->skill_lv + m->skill_base) / 100.0;
    if(m->skill_target == self_buff && 
        (m->skill_buff == buff_atk || m->skill_buff == buff_atk_def))
      return m->attack * str;
    else {
      if(m->skill_target == team_buff &&
          (m->skill_buff == buff_atk || m->skill_buff == buff_atk_def)){
        double r = 0.0;
        if(m->skill_attr.empty()){
          for(const auto& c: our){ r += c->attack * str; }
        } else {
          for(const auto& c: our){
            if(c->attr == m->skill_attr)
              r += c->attack * str;
          }
        }
        return r;
      }
      else if(m->skill_target == enemy_buff &&
          (m->skill_buff == buff_def || m->skill_buff == buff_atk_def)){
        if(m->skill_attr.empty()) return defense * str;
        else if(m->skill_attr == attr) return defense * str;
        else return 0.0;
      }
      else return 0.0;
    }
  }
}

std::vector<std::pair<double, uint64_t>> compute_skill_distribution(
    const std::vector<std::pair<double, int>> skills){
  std::vector<std::pair<double, uint64_t>> res;
  for(const auto& s: selections){
    double skill = 0.0;
    uint64_t prob = 1;
    for(int i = 0; i < 5; ++i){
      if(s[i]){
        skill += skills[i].first;
        prob *= skills[i].second;
      } else prob *= (100 - skills[i].second);
    }
    res.push_back(std::make_pair(skill, prob));
  }
  return res;
}

template <typename E>
int try_skill_attack(const std::vector<std::shared_ptr<monster>>& our,
    const std::shared_ptr<E> target){
  int base = 0;
  for(const auto& m : our) { base += m->attack; }

  double skill_atk = 0.0;
  for(const auto& m : our) { skill_atk += compute_skill(our, m, target->defense, target->attr); }

  for(int i = 0; i < 5; i++){
    auto atk = compute_basic_atk(our, target, base, i) + skill_atk;
    if(atk > target->defense) return i;
  }
  return -1;
}

template <typename E>
std::pair<int, double> highest_chance(const std::vector<std::shared_ptr<monster>>& our,
    std::shared_ptr<E> enemy){
  auto chance = 0.0;
  int base = 0;
  int ld = 0;
  for(const auto& m : our) { base += m->attack; }

  for(int i = 0; i < 5; ++i){
    auto leader = our[i];
    std::vector<std::pair<double, int>> skills;
    for(int j = 0; j < 5; ++j){
      const auto& m = our[j];
      auto sk = compute_skill(our, m, enemy->defense, enemy->attr);
      auto freq = (i == j) ? m->skill_freq + 15 : m->skill_freq;
      skills.push_back(std::make_pair(sk, freq));
    }

    auto dist = compute_skill_distribution(skills);
    auto threshold = enemy->defense - compute_basic_atk(our, enemy, base, i);

    uint64_t win_chance = 0;
    for(const auto& p: dist){
      if(p.first > threshold)  win_chance += p.second;
    }

    double f_chance = win_chance / 10000000000.0;
    if(f_chance > chance) {
      chance = f_chance;
      ld = i;
    }
  }

  return std::make_pair(ld, chance);
}

template <class E>
class nonskill_f{
  const std::vector<std::shared_ptr<E>>& enemys;

  public:
  std::vector<double> chance;

  std::vector<std::vector<std::shared_ptr<monster>>> result;

  nonskill_f(const std::vector<std::shared_ptr<E>>& e): enemys(e), result(e.size()), chance(e.size(), 0.0){};

  template <class It>
    bool operator()(It first, It last){
      auto ms = std::vector<std::shared_ptr<monster>>(first, last);
      int cost = 0;
      for(const auto& m : ms) cost += m->cost;

      auto r_it = result.begin();
      for(auto e_it = enemys.begin(), e_end = enemys.end(); 
          e_it != e_end; ++e_it, ++r_it){
        if(try_attack(ms, *e_it) != -1){
          if(r_it->size() == 0) *r_it = ms;
          else {
            int o_cost = 0;
            for(const auto& m : (*r_it)) o_cost += m->cost;
            if(o_cost > cost) *r_it = ms;
          }
        }
      }
      return false;
    }
};

template<typename S, typename E>
void printTeam(S& os, const std::vector<std::shared_ptr<monster>>& team, double atk, E e){
  int cost = 0;
  for(const auto& m : team) cost += m->cost;
  os << e->name << std::endl;
  os << team[0]->id;
  for(auto i = team.begin() + 1, e = team.end(); i != e; ++i){
    os << ',' << (*i)->id;
  }
  os << std::endl;
  os << "cost:" << cost << std::endl;
  os << team[0]->name << '(' << team[0]->lv << ')';
  for(auto i = team.begin() + 1, e = team.end(); i != e; ++i){
    os << ',' << (*i)->name << '(' << (*i)->lv << ')';
  }
  os << std::endl;
  os << team[0]->attack;
  for(auto i = team.begin() + 1, e = team.end(); i != e; ++i){
    os << ',' << (*i)->attack;
  }
  os << " => " << atk << " > " << e->defense;
}

template <class E>
class skill_f{
  const std::vector<std::shared_ptr<E>>& enemys;

  public:
  std::vector<std::vector<std::shared_ptr<monster>>> result;

  std::vector<double> chance;

  skill_f(const std::vector<std::shared_ptr<E>>& e): enemys(e), result(e.size()), chance(e.size(), 0.0){};

  template <class It>
    bool operator()(It first, It last){
      auto ms = std::vector<std::shared_ptr<monster>>(first, last);
      int cost = 0;
      for(const auto& m : ms) cost += m->cost;

      auto r_it = result.begin();
      auto c_it = chance.begin();
      for(auto e_it = enemys.begin(), e_end = enemys.end(); 
          e_it != e_end; ++e_it, ++r_it, ++c_it){
        if(try_attack(ms, *e_it) != -1){
          if(r_it->size() == 0) {
            auto chance = highest_chance(ms, *e_it);
            *r_it = ms;
            *c_it = chance.second;
          }
          else {
            int o_cost = 0;
            for(const auto& m : (*r_it)) o_cost += m->cost;

            auto chance = highest_chance(ms, *e_it);

            if((chance.second > *c_it) || (chance.second == *c_it && o_cost > cost)){
              *r_it = ms;
              *c_it = chance.second;
            }
          }
        }
      }
      return false;
    }
};

template <typename E>
void nonskill(std::vector<std::shared_ptr<monster>> our, 
    const std::vector<std::shared_ptr<E>>& enemys,
    const std::string& prefix, const std::string& postfix){
  std::ofstream ofs(prefix + "_" + postfix + ".txt", std::ofstream::trunc);
  if(!ofs.is_open()){
      std::cerr << "Error opening file: " << prefix + "_" + postfix + ".txt" << std::endl;
      return;
  }
  auto r = for_each_combination(our.begin(), our.begin() + 5, our.end(), nonskill_f<E>(enemys));
  auto result = r.result;

  for(int i = 0, j = result.size(); i < j; ++i){
    auto e = enemys[i];
    auto r = result[i];
    if(r.size() != 0){
      auto idx = try_attack(r, e);
      auto leader = r[idx];
      std::vector<std::shared_ptr<monster>> team;
      team.push_back(leader);
      team.insert(team.end(), r.begin(), r.begin() + idx);
      if(idx != 4)
        team.insert(team.end(), r.begin() + idx + 1, r.end());

      int base = 0;
      for(const auto& m : team) base += m->attack; 

      printTeam(ofs, team, compute_basic_atk(team, e, base, 0), e);
      ofs << std::endl << std::endl;

    } else {
      ofs << e->name << ": no team" << std::endl << std::endl;
    }
  }

  ofs.close();
  return;
}

template <typename E>
void skill(std::vector<std::shared_ptr<monster>> our, 
    const std::vector<std::shared_ptr<E>>& enemys,
    const std::string& prefix, const std::string& postfix){
  std::ofstream ofs(prefix + "_" + postfix + ".txt", std::ofstream::trunc);
  if(!ofs.is_open()){
      std::cerr << "Error opening file: " << prefix + "_" + postfix + ".txt" << std::endl;
      return;
  }
  auto r = for_each_combination(our.begin(), our.begin() + 5, our.end(), skill_f<E>(enemys));
  auto result = r.result;

  for(int i = 0, j = result.size(); i < j; ++i){
    auto e = enemys[i];
    auto r = result[i];
    if(r.size() != 0){
      auto idx_ = highest_chance(r, e);
      auto idx = idx_.first;
      auto leader = r[idx];
      std::vector<std::shared_ptr<monster>> team;
      team.push_back(leader);
      team.insert(team.end(), r.begin(), r.begin() + idx);
      if(idx != 4)
        team.insert(team.end(), r.begin() + idx + 1, r.end());

      int base = 0;
      for(const auto& m : team) base += m->attack; 

      double skill_atk = 0.0;
      for(const auto& m : r) { skill_atk += compute_skill(r, m, e->defense, e->attr); }

      printTeam(ofs, team, compute_basic_atk(team, e, base, 0) + skill_atk, e);
      ofs << std::endl;
      ofs << "prob: " << idx_.second << std::endl;
      ofs << std::endl << std::endl;

    } else {
      ofs << e->name << ": no team" << std::endl << std::endl;
    }
  }

  ofs.close();
  return;
}

template <typename E>
int parse_file(const std::string& file,
    std::shared_ptr<E> (*line_parser)(const std::string&),
    std::vector<std::shared_ptr<E>>& out){
    std::ifstream fmonster;
    fmonster.open(file);
    if(fmonster.is_open()){
      int l = 1;
      while(!fmonster.eof()){
        std::string line;
        std::getline(fmonster, line);
        auto m = line_parser(line);
        if(m){
          out.push_back(m);
        } else if(!line.empty()){
          std::cerr << "Parsing error: Line" << l << " of " << file << std::endl;
          return EXIT_FAILURE;
        }
        ++l;
      }
    } else {
      std::cerr << "Error opening file " << file << std::endl;
      return EXIT_FAILURE;
    }
    fmonster.close();
    std::cout << out.size() << " data loaded." << std::endl;
    return EXIT_SUCCESS;
}

int main (int argc, char *argv[]){
  if(argc < 2){
      std::cerr << "Please input file prefix." << std::endl;
      return EXIT_FAILURE;
  }

  std::string prefix;
  if(argc >= 2) prefix.assign(argv[1]);

  bool enable[] = {true, true, true, true};
  if(argc > 2){
      enable[0] = enable[1] = enable[2] = enable[3] = false;
  }

  const std::string& enable_opt("enable");
  for(int i = 0; i < (argc - 2) && i < 4; ++i){
      if(i < (argc - 2) && enable_opt == argv[i + 2])
          enable[i] = true;
  }

  std::vector<std::shared_ptr<monster>> monster;
  std::vector<std::shared_ptr<enemy>> basic_enemys;
  std::vector<std::shared_ptr<enemy_reduce_attack>> reduce_enemys;
  std::vector<std::shared_ptr<enemy>> skill_enemys;
  std::vector<std::shared_ptr<enemy_reduce_attack>> reduce_skill_enemys;

  std::cout << "Loading daemon data." << std::endl;
  if(parse_file(prefix + "_daemons.csv", *parse_monster, monster)
      == EXIT_FAILURE)
  { return EXIT_FAILURE; }

  std::cout << "Loading basic monsters data." << std::endl;
  if(parse_file("tower_basic.csv", *parse_enemy, basic_enemys)
      == EXIT_FAILURE)
  { basic_enemys.clear(); }

  std::cout << "Loading data of monsters which need skill." << std::endl;
  if(parse_file("tower_need_skill.csv", *parse_enemy, skill_enemys)
      == EXIT_FAILURE)
  { skill_enemys.clear(); }

  std::cout << "Loading data of monsters which reduce our attack." << std::endl;
  if(parse_file("tower_reduce_atk.csv", *parse_enemy_reduce_attack, reduce_enemys)
      == EXIT_FAILURE)
  { reduce_enemys.clear(); }

  std::cout << "Loading data of monsters which reduce our attack and need skill." << std::endl;
  if(parse_file("tower_reduce_need_skill.csv", *parse_enemy_reduce_attack, reduce_skill_enemys)
      == EXIT_FAILURE)
  { reduce_skill_enemys.clear(); }

  auto start = std::chrono::system_clock::now();

  if(enable[0] && basic_enemys.size() != 0) 
      nonskill(monster, basic_enemys, prefix, "tower_guide"); 
  else std::cout << "Skipping basic monsters." << std::endl;

  if(enable[1] && reduce_enemys.size() != 0) 
  nonskill(monster, reduce_enemys, prefix, "tower_guide_reduce");
  else std::cout << "Skipping monsters which reduce team's attack." << std::endl;

  if(enable[2] && skill_enemys.size() != 0) 
  skill(monster, skill_enemys, prefix, "tower_guide_need_skill");
  else std::cout << "Skipping monsters which need daemon's skill." << std::endl;

  if(enable[3] && reduce_skill_enemys.size() != 0) 
  skill(monster, reduce_skill_enemys, prefix, "tower_guide_reduce_need_skill");
  else std::cout << "Skipping monsters which reduce team's attack and needs daemon's skill." << std::endl;
  auto end = std::chrono::system_clock::now();
  std::cout << "Took "
   << (end - start) / std::chrono::seconds(1)
   << "seconds." << std::endl;
  return EXIT_SUCCESS;
}

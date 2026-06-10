// Doppio-port phase 3 bridge (2026-06-06): byte-identical CLI vs C ABI
// snapshot for find_lr_orders_scan.
//
// Locks the transport-neutral contract of the projective Cheng-Wu
// gauge-scan op: both `hyperflint eval-json` (CLI) and
// `hf_find_lr_orders_scan` (C ABI) route through
// hyperflint::handlers::find_lr_orders_scan, which emits the
// {schema_version, hf_version} envelope inline (so splice_envelope
// returns the body unchanged — branch 1, same as find_lr_orders).
//
// Cases:
//   A* — ABI sanity (non-NULL, single envelope)
//   B* — CLI sanity (non-empty, op-prefixed, envelope present)
//   C1 — byte-identity modulo timing_compute_s
//   D1 — missing-exps error path identical on both transports
//   E1 — deterministic scan content on the massless-box fixture
//        (projective, >= 1 order, first gauge/order stable)
//
// Mirrors test_c_abi_find_lr_orders_cli_snapshot.cpp (iter-58); the
// strip_timing / run_cli helpers are intentionally duplicated, per the
// "snapshot tests are independent" lineage expectation.
//
// Usage:
//   test_c_abi_find_lr_orders_scan_cli_snapshot <path-to-hyperflint-cli>

#include "hyperflint/c_abi.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>  // close, unlink, mkstemp

namespace {

std::string strip_timing(const std::string& s) {
    const std::string marker = ",\"timing_compute_s\":";
    size_t a = s.find(marker);
    if (a == std::string::npos) return s;
    size_t b = a + marker.size();
    auto is_num_char = [](char c) {
        return (c >= '0' && c <= '9') || c == '.' || c == '-' || c == '+'
            || c == 'e' || c == 'E'
            || c == 'n' || c == 'N' || c == 'a' || c == 'A'
            || c == 'i' || c == 'I' || c == 'f' || c == 'F'
            || c == '(' || c == ')' || c == '_';
    };
    while (b < s.size() && is_num_char(s[b])) ++b;
    return s.substr(0, a) + s.substr(b);
}

std::string run_cli(const std::string& cli_path,
                    const std::string& request) {
    char req_path[] = "/tmp/hf_cli_lrscan_snapshot_req_XXXXXX";
    char out_path[] = "/tmp/hf_cli_lrscan_snapshot_out_XXXXXX";
    int rfd = ::mkstemp(req_path);
    int ofd = ::mkstemp(out_path);
    if (rfd < 0 || ofd < 0) {
        std::cerr << "[FAIL] mkstemp failed\n";
        return "";
    }
    ::close(rfd);
    ::close(ofd);
    {
        std::ofstream f(req_path);
        f << request;
    }
    std::string cmd = cli_path + " eval-json < " + req_path
                    + " > " + out_path + " 2>/dev/null";
    int rc = std::system(cmd.c_str());
    if (rc != 0) {
        std::cerr << "[WARN] CLI returned exit code " << rc << "\n";
    }
    std::ifstream f(out_path);
    std::stringstream ss;
    ss << f.rdbuf();
    ::unlink(req_path);
    ::unlink(out_path);
    std::string out = ss.str();
    if (!out.empty() && out.back() == '\n') out.pop_back();
    return out;
}

size_t count_substr(const std::string& hay, const std::string& need) {
    size_t n = 0, pos = 0;
    while ((pos = hay.find(need, pos)) != std::string::npos) {
        ++n;
        pos += need.size();
    }
    return n;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <path-to-hyperflint-cli>\n";
        return 2;
    }
    const std::string cli_path = argv[1];

    int n_pass = 0;
    int n_fail = 0;
    auto check = [&](bool cond, const char* tag) {
        if (cond) { std::printf("[PASS] %s\n", tag); ++n_pass; }
        else      { std::printf("[FAIL] %s\n", tag); ++n_fail; }
    };

    // Massless box U+F scan fixture: projective at exps
    // U^(0+2eps) F^(-2-eps): sum a_i d_i = 0*1 + (-2)*2 = -4 = -n,
    // sum b_i d_i = 2*1 + (-1)*2 = 0.  Sub-second; deterministic.
    const std::string req =
        "{\"op\":\"find_lr_orders_scan\","
         "\"groups\":[[\"x1 + x2 + x3 + x4\",\"s*x1*x3 + t*x2*x4\"]],"
         "\"xvars\":[\"x1\",\"x2\",\"x3\",\"x4\"],"
         "\"coeff_vars\":[\"s\",\"t\"],"
         "\"exps\":[[[0,2],[-2,-1]]],"
         "\"keep_rule\":\"Strict\","
         "\"max_orders\":64}";

    // ---- C ABI ----
    char* abi_raw = hf_find_lr_orders_scan(req.c_str());
    check(abi_raw != nullptr, "A1 hf_find_lr_orders_scan: non-NULL return");
    std::string abi_full = abi_raw ? std::string(abi_raw) : std::string();
    hf_free_string(abi_raw);
    check(count_substr(abi_full, "\"schema_version\":") == 1,
          "A2 ABI: schema_version exactly once (no double-wrap)");
    check(count_substr(abi_full, "\"hf_version\":") == 1,
          "A3 ABI: hf_version exactly once (no double-wrap)");
    std::string abi_stripped = strip_timing(abi_full);

    // ---- CLI ----
    std::string cli_out = run_cli(cli_path, req);
    check(!cli_out.empty(), "B1 CLI: non-empty stdout");
    check(cli_out.rfind("{\"op\":\"find_lr_orders_scan\"", 0) == 0,
          "B2 CLI: output starts with {\"op\":\"find_lr_orders_scan\"");
    check(cli_out.find("\"schema_version\":") != std::string::npos,
          "B3 CLI: envelope present (handler emits inline)");
    std::string cli_stripped = strip_timing(cli_out);

    // ---- byte identity ----
    bool equal = (abi_stripped == cli_stripped);
    check(equal, "C1 byte-identical CLI vs C-ABI (timing_compute_s stripped)");
    if (!equal) {
        std::cerr << "    CLI stripped: " << cli_stripped << "\n";
        std::cerr << "    ABI stripped: " << abi_stripped << "\n";
    }

    // ---- missing-exps error path, identical on both transports ----
    const std::string bad_req =
        "{\"op\":\"find_lr_orders_scan\","
         "\"groups\":[[\"x1 + x2\",\"s*x1*x2\"]],"
         "\"xvars\":[\"x1\",\"x2\"],\"coeff_vars\":[\"s\"]}";
    char* bad_abi_raw = hf_find_lr_orders_scan(bad_req.c_str());
    std::string bad_abi = bad_abi_raw ? std::string(bad_abi_raw)
                                      : std::string();
    hf_free_string(bad_abi_raw);
    std::string bad_cli = run_cli(cli_path, bad_req);
    check(bad_abi.find("\"error\":") != std::string::npos,
          "D1a ABI: missing exps -> error envelope");
    check(strip_timing(bad_abi) == strip_timing(bad_cli),
          "D1b missing-exps error byte-identical CLI vs ABI");

    // ---- deterministic scan content ----
    check(abi_full.find("\"projective\":true") != std::string::npos,
          "E1a box fixture: projective");
    check(abi_full.find("\"orders\":[{") != std::string::npos,
          "E1b box fixture: at least one admissible (gauge, order)");
    check(abi_full.find("\"carried_sqrts\":0") != std::string::npos,
          "E1c box fixture: strict rule => no carried sqrts");

    std::printf("[SUMMARY] PASS=%d FAIL=%d\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}

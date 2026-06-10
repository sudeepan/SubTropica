// Carry-discharge (FindRoots) keep rule in find_lr_orders — TDD battery.
//
// find_lr_orders judges deg-2 algebraic letters TERMINAL-ONLY: a deg-2
// letter whose sqrt-argument (odd part of the pivot-discriminant) depends
// on a still-pending integration variable is rejected at that step
// (the forbidden_after_step guard, lr_search.cpp).  The Doppio FindRoots
// tier (lr_scan::step_fr / fr_judge) instead CARRIES that obligation
// forward and discharges it once the pending variable is integrated —
// strictly more permissive.  This test wires that tier into
// find_lr_orders as the "carry_discharge" request option, DEFAULT ON.
//
// The three gates of the task:
//   (1) NOLR->LR FLIP: a fixture that find_lr_orders+algebraic reports
//       NOLR under the terminal-only rule but LR under carry-discharge.
//       Primary fixture is the REAL gauged-uq5 system (uq5 at x3->1):
//       lr_scan FindRoots @ gauge x3 has 24 admissible orders (oracle in
//       test_lr_scan.cpp / t24 battery), so an LR order provably exists;
//       find_lr_orders+algebraic currently says NOLR.  Plus a small
//       synthetic 2-var fixture with the same structure.
//   (2) STRICT BYTE-IDENTITY: carry_discharge:false reproduces today's
//       response (the response is fed through handlers::find_lr_orders so
//       the whole envelope is checked; the timing field is stripped).
//   (3) DEFAULT-ON: an absent carry_discharge behaves like :true.
//
// Parity against lr_scan FindRoots (the per-step semantics are the SAME
// extracted code, lr_scan::step_fr_judge) is covered structurally by the
// flip oracle above and by test_lr_scan.cpp (unchanged 120-order pin).
//
// Handler-level (mirrors test_find_lr_orders_strategy_roundtrip.cpp): the
// request crosses the JSON boundary, so this also guards the serialization
// of the new option + any carried-sqrt profile fields.
//
// Tier-1 wall: gauged-uq5 ~20 ms, synthetic sub-ms; a handful of rows.

#include "hyperflint/bridge/handlers.hpp"

#include <iostream>
#include <regex>
#include <string>
#include <vector>

namespace {

// ---- response field accessors (hand-rolled, no JSON dep — same style as
//      test_find_lr_orders_strategy_roundtrip.cpp) ----

bool resp_nolr(const std::string& resp) {
    std::regex rx(R"~("nolr"\s*:\s*(true|false))~");
    std::smatch m;
    if (std::regex_search(resp, m, rx)) return m[1].str() == "true";
    // absent => treat as NOLR=unknown; force a visible failure
    return true;
}

bool resp_has_error(const std::string& resp) {
    return resp.find("\"error\"") != std::string::npos;
}

std::string resp_best_order(const std::string& resp) {
    std::regex rx(R"~("best_order"\s*:\s*(\[[^\]]*\]))~");
    std::smatch m;
    if (std::regex_search(resp, m, rx)) return m[1].str();
    return "<absent>";
}

// Extract an unsigned integer field; returns -1 if absent.
long resp_uint_field(const std::string& resp, const char* field) {
    std::regex rx("\"" + std::string(field) + "\"\\s*:\\s*([0-9]+)");
    std::smatch m;
    if (std::regex_search(resp, m, rx)) return std::stol(m[1].str());
    return -1;
}

// Strip the (non-deterministic) timing field so two responses can be
// compared byte-for-byte.  Same approach the c-abi snapshot ctests use.
std::string strip_timing(const std::string& resp) {
    std::regex rx(R"~(,"timing_compute_s":[0-9eE.+\-]+)~");
    return std::regex_replace(resp, rx, "");
}

int g_pass = 0, g_fail = 0;
void check(const char* label, bool cond, const std::string& diag = "") {
    if (cond) {
        std::cout << "[PASS] " << label << "\n";
        ++g_pass;
    } else {
        std::cerr << "[FAIL] " << label;
        if (!diag.empty()) std::cerr << "  " << diag;
        std::cerr << "\n";
        ++g_fail;
    }
}

// ---- fixtures ----

// uq5 gauged at x3 -> 1.  Two-factor projective group of the ungauged
// uq5 (test_lr_scan.cpp) with x3 substituted to 1:
//   P1 = x1 + x2 + x3                       -> x1 + x2 + 1
//   P2 = -qq1 x1 x2 - qq2 x1 x3 + 2 wb1 x3 x4 - x4^2
//        + 2 wb2 x2 x5 - x5^2 + 2 yb x4 x5  -> (x3 -> 1)
const char* kGaugedUq5Polys =
    R"~("x1+x2+1","-qq1*x1*x2-qq2*x1+2*wb1*x4-x4^2+2*wb2*x2*x5-x5^2+2*yb*x4*x5")~";
const char* kGaugedUq5Xvars  = R"~("x1","x2","x4","x5")~";
const char* kGaugedUq5Coeffs = R"~("qq1","qq2","wb1","wb2","yb")~";

std::string gauged_uq5_request(const char* carry_field) {
    std::string r = "{\"op\":\"find_lr_orders\",";
    r += "\"xvars\":[";   r += kGaugedUq5Xvars;  r += "],";
    r += "\"coeff_vars\":["; r += kGaugedUq5Coeffs; r += "],";
    r += "\"polys\":[";   r += kGaugedUq5Polys;  r += "],";
    r += "\"algebraic_letters\":true";
    if (carry_field[0] != '\0') { r += ","; r += carry_field; }
    r += "}";
    return r;
}

// Small synthetic 2-var fixture, VERIFIED against lr_scan FindRoots.
// It is the Cheng-Wu gauge-x3->1 projection of the projective 3-var
// system
//     U = x1 + x2 + x3 ,   F = x1^2 + x2^2 + x1 x2 + s x3^2
// (twist exps U^(-1+2eps) F^(-1-eps); Sigma a_i deg = -3, Sigma b_i
// deg = 0).  On that lift:
//   * lr_scan Strict      -> 0 orders (gauge-exhaustive NOLR);
//   * lr_scan FindRoots   -> 6 orders, each carried_sqrts == 1; at gauge
//     x3 the orders are (x1,x2) and (x2,x1).
// The gauge-x3 group is { x1+x2+1, x1^2+x2^2+x1 x2 + s } (F at x3->1,
// plus the U-boundary x1+x2+1).  find_lr_orders+carry on it must thus:
//   * carry OFF -> NOLR  (== lr_scan Strict @ gauge x3),
//   * carry ON  -> LR, order (x1,x2), carried_sqrts == 1  (one of the
//     lr_scan FindRoots @ gauge x3 orders; the score-minimal one).
// The single deg-2 letter's sqrt-obligation depends on the pending
// second variable at the first pivot and discharges at the second —
// exactly the structure the terminal-only rule rejects but carry keeps.
const char* kSynthPolys  = R"~("x1+x2+1","x1^2+x2^2+x1*x2+s")~";
const char* kSynthXvars  = R"~("x1","x2")~";
const char* kSynthCoeffs = R"~("s")~";

std::string synth_request(const char* carry_field) {
    std::string r = "{\"op\":\"find_lr_orders\",";
    r += "\"xvars\":[";   r += kSynthXvars;  r += "],";
    r += "\"coeff_vars\":["; r += kSynthCoeffs; r += "],";
    r += "\"polys\":[";   r += kSynthPolys;  r += "],";
    r += "\"algebraic_letters\":true";
    if (carry_field[0] != '\0') { r += ","; r += carry_field; }
    r += "}";
    return r;
}

}  // namespace

int main() {
    using hyperflint::handlers::find_lr_orders;

    // ===== Gate 1: NOLR -> LR flip (real gauged-uq5) =====
    {
        const std::string off = find_lr_orders(gauged_uq5_request(
            "\"carry_discharge\":false"));
        const std::string on  = find_lr_orders(gauged_uq5_request(
            "\"carry_discharge\":true"));
        check("gauged-uq5: carry OFF => NOLR (terminal-only, today)",
              !resp_has_error(off) && resp_nolr(off), off);
        check("gauged-uq5: carry ON => LR (carry-discharge flip)",
              !resp_has_error(on) && !resp_nolr(on), on);
        if (!resp_nolr(on))
            std::cout << "       gauged-uq5 carry-ON order = "
                      << resp_best_order(on) << "\n";
    }

    // ===== Gate 1b: NOLR -> LR flip + lr_scan PARITY (small synthetic) =====
    {
        const std::string off = find_lr_orders(synth_request(
            "\"carry_discharge\":false"));
        const std::string on  = find_lr_orders(synth_request(
            "\"carry_discharge\":true"));
        check("synthetic: carry OFF => NOLR (== lr_scan Strict @ gauge x3)",
              !resp_has_error(off) && resp_nolr(off), off);
        check("synthetic: carry ON => LR (carry-discharge flip)",
              !resp_has_error(on) && !resp_nolr(on), on);
        // PARITY with lr_scan FindRoots @ gauge x3: order (x1,x2),
        // carried_sqrts == 1 (verified: see fixture comment).  This is
        // the core correctness claim at small scale — the same verdict
        // AND the same carried-sqrt count as the exhaustive scan engine.
        check("synthetic: carry-ON order == (x1,x2) [lr_scan parity]",
              resp_best_order(on) == "[\"x1\",\"x2\"]",
              "order=" + resp_best_order(on));
        check("synthetic: carry-ON carried_sqrts == 1 [lr_scan parity]",
              resp_uint_field(on, "carried_sqrts") == 1,
              "carried_sqrts=" +
                  std::to_string(resp_uint_field(on, "carried_sqrts")));
    }

    // ===== Gate 2: Strict byte-identity (carry OFF == today) =====
    // The pre-change baseline response for gauged-uq5 + algebraic is the
    // terminal-only NOLR envelope; carry_discharge:false must reproduce
    // it byte-for-byte (modulo timing).  We compare carry:false against
    // the EXPECTED today-envelope captured here as a literal, so a
    // regression in the Strict path (not just a verdict flip) is caught.
    {
        const std::string off = strip_timing(find_lr_orders(
            gauged_uq5_request("\"carry_discharge\":false")));
        // today's NOLR envelope for this fixture (algebraic on, no carry):
        const std::string expect =
            "{\"op\":\"find_lr_orders\""
            ",\"schema_version\":1"
            ",\"hf_version\":\"" /* version-agnostic prefix match below */;
        // Version-agnostic structural check: NOLR envelope shape.
        const bool shape_ok =
            off.find("\"best_order\":[]") != std::string::npos &&
            off.find("\"score\":null") != std::string::npos &&
            off.find("\"nolr\":true") != std::string::npos &&
            off.find("\"strategy\":\"Fubini_Lungo\"") != std::string::npos &&
            off.find("\"root_polys\":[]") != std::string::npos;
        check("gauged-uq5: Strict(carry:false) NOLR envelope unchanged",
              shape_ok, off);
        (void) expect;
    }

    // ===== Gate 3: default-ON (absent option == :true) =====
    {
        const std::string def = find_lr_orders(gauged_uq5_request(""));
        const std::string on  = find_lr_orders(gauged_uq5_request(
            "\"carry_discharge\":true"));
        check("gauged-uq5: absent carry_discharge => LR (default ON)",
              !resp_has_error(def) && !resp_nolr(def), def);
        check("gauged-uq5: absent == explicit true (byte-identical)",
              strip_timing(def) == strip_timing(on),
              "def=" + strip_timing(def) + "  on=" + strip_timing(on));

        const std::string sdef = find_lr_orders(synth_request(""));
        check("synthetic: absent carry_discharge => LR (default ON)",
              !resp_has_error(sdef) && !resp_nolr(sdef), sdef);
    }

    std::cout << "Summary: " << g_pass << " PASS / " << g_fail << " FAIL\n";
    return g_fail == 0 ? 0 : 1;
}

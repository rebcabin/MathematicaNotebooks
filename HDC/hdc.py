r"""By Brian Beckman and Dan Ilan

MIT licensed -- see LICENSE.txt in this repository.

https://www.hd-computing.com/

https://www.dropbox.com/scl/fi/2nqq8qhmr8y1l1jlvsllz/JoshiEtAl-QI2016-language-geometry-copy.pdf?rlkey=sx1htkhikfar2cary8pqq3yhu&dl=0

https://bitbucket.org/gsitech/language_indexing/src/master/?search_id=4f2d14b1-c9af-486f-9dc4-4333f5fcae37

Algorithm

The GSI Hyperdimensional Computing (HDC) program is used to
determine the language of text using the algorithm described
in Language Geometry Hypervector Embedding.

To do so, the program:

    Trains binary Hyper-Dimensional (HD) vectors for each
    language from a set of texts (training set) in that
    language.

    Every character in the English alphabet and the space
    character is assigned a randomly generated HD binary
    vector.

    During training, N-grams are generated from texts (see
    paper for details), converted to bi-polar (‘1’ or ‘-1’)
    and accumulated in a signed HD vector.

    After processing all N-grams of a text, or set of texts,
    the accumulator vector is converted to a binary vector by
    taking the sign of each element in the vector. The result
    is the reference vector of the text or language.

    Inference, or classification, is done by computing the
    Hamming distance between the trained vector of a query
    text, to all vectors of the trained languages. The
    language that produces the minimal Hamming distance is
    predicted to be the query text language. Note: This step
    varies from the process described in the paper.

The document has pointers to the code in our bitbucket git
repository.

"""

from belex.literal import (GL, NOOP, RL, VR, WRL, ERL, RSP16,
                           Mask, Section, apl_commands,
                           belex_apl, u16)

from belex_libs.common import cpy_imm_16
from belex_libs.arithmetic import add_u16

from belex.utils.example_utils import u16_to_vr


#   _   _ _   _ _ _ _   _
#  | | | | |_(_) (_) |_(_)___ ___
#  | |_| |  _| | | |  _| / -_|_-<
#   \___/ \__|_|_|_|\__|_\___/__/


@belex_apl
def hdc_cp(Belex,
           dst: VR, d: Section,
           src: VR, s: Section) -> None:
    RL[s]  <= src()
    GL[s]  <= RL()
    dst[d] <= GL()


@belex_apl
def hdc_swap(Belex,
           dst: VR, a: Section,
           src: VR, b: Section) -> VR:
    """a = b ^ a; b = b ^ a; a = b ^ a"""
    # TODO
    pass


#   __  __      _            __  _____  ___
#  |  \/  |_  _| |  ___ ___  \ \/ / _ \| _ \
#  | |\/| | || | | |___|___|  >  < (_) |   /
#  |_|  |_|\_,_|_| |___|___| /_/\_\___/|_|_\


@belex_apl
def hdc_mul_const_vr (Belex,
                      dst : VR,
                      d   : Section,
                      src : VR,
                      x   : Section,
                      y   : Section)  -> None:
    """In hdc, mul is xor."""
    RL[x]  <= src()
    GL[x]  <= RL()
    RL[y]  <= src() ^ GL()
    GL[y]  <= RL()
    dst[d] <= GL()


@belex_apl
def hdc_mul_const_section (Belex,
                           dst : VR,
                           x   : VR,
                           y   : VR,
                           sec : Section) -> None:
    """In hdc, mul is xor."""
    RL[sec]  <= x()
    RL[sec]  ^= y()
    dst[sec] <= RL()


#   ____    __  __       _         _ _
#  |__ /___|  \/  |__ _ (_)___ _ _(_) |_ _  _
#   |_ \___| |\/| / _` || / _ \ '_| |  _| || |
#  |___/   |_|  |_\__,_|/ \___/_| |_|\__|\_, |
#                     |__/               |__/


@belex_apl
def hdc_ternary_majority_const_section (
        Belex,
        dst : VR,
        a   : VR,
        b   : VR,
        c   : VR,
        sec : Section) -> None:
    """In hdc, ADD is majority. This routine works
    only if all VRs are different.
    """
    RL[sec]  <= a() & b()
    RL[sec]  <= b() & c() | RL()
    RL[sec]  <= c() & a() | RL()
    dst[sec] <= RL()


@belex_apl
def hdc_ternary_majority_groundpound_all_sections(
        Belex,
        dst : VR,
        a   : VR,
        b   : VR,
        c   : VR,
) -> None:

    def groundpound(sec : str):
        RL[sec]  <= a() & b()
        RL[sec]  <= b() & c() | RL()
        RL[sec]  <= c() & a() | RL()
        dst[sec] <= RL()

    for sec in ["F", "E", "D", "C", "B", "A", "9", "8",
                "7", "6", "5", "4", "3", "2", "1", "0"]:
        groundpound(sec)
    pass


@belex_apl
def hdc_ternary_majority_const_vr (
        Belex,
        dst : VR, d : Section,
        src : VR, a : Section,
                  b : Section,
                  c : Section) -> None:
    """Compute majority for three sections in a
    single VR and propagate it to a section d in
    dst, a different VR. d may not equal any of the
    other sections as src[d] is scratch space. That
    fits with our assumption that at most 15
    sections of a single VR are live, with one
    reserved for scratch.
    """
    RL[a] <= src()
    GL[a] <= RL()
    RL[b] <= src() & GL()
    GL[b] <= RL()
    RL[d] <= GL()

    RL[b] <= src()
    GL[b] <= RL()
    RL[c] <= src() & GL()
    GL[c] <= RL()
    RL[d] |= GL()

    RL[c] <= src()
    GL[c] <= RL()
    RL[a] <= src() & GL()
    GL[a] <= RL()
    RL[d] |= GL()

    dst[d] <= RL()
    pass


#   ___     __  __       _         _ _
#  / _ \___|  \/  |__ _ (_)___ _ _(_) |_ _  _
#  \_, /___| |\/| / _` || / _ \ '_| |  _| || |
#   /_/    |_|  |_\__,_|/ \___/_| |_|\__|\_, |
#                     |__/               |__/


@belex_apl
def hdc_9_maj_via_3_const_section(
        Belex,
        dv : VR,
        a : VR, b : VR, c : VR,
        d : VR, e : VR, f : VR,
        g : VR, h : VR, i : VR,
        s : Section) -> None:
    """27 fewer instructions than hdc_9_maj_via_3_const_vr, plus
    does not require two fragments. Uses sections 0, 1, 2, 3, 4, 5
    of tmp1, so the input section, s, cannot be any of these. It
    is possible to engineer around this limitation, but might not
    be worth it at the moment."""
    tmp1 : VR = Belex.VR()
    tmp2 : VR = Belex.VR()

    t1mabc : Section = Belex.Section(0)
    hdc_ternary_majority_const_section(tmp1, a, b, c, s)
    hdc_cp(tmp1, t1mabc, tmp1, s)

    t1mdef : Section = Belex.Section(1)
    hdc_ternary_majority_const_section(tmp1, d, e, f, s)
    hdc_cp(tmp1, t1mdef, tmp1, s)

    t1mghi : Section = Belex.Section(2)
    hdc_ternary_majority_const_section(tmp1, g, h, i, s)
    hdc_cp(tmp1, t1mghi, tmp1, s)

    hdc_cp(tmp2, "0", tmp1, t1mdef)
    hdc_cp(tmp2, "1", c, s)
    hdc_cp(tmp2, "2", b, s)
    hdc_cp(tmp2, "3", a, s)

    # l__("4") = 3maj(c, mdef, b)
    hdc_ternary_majority_const_vr(tmp2, "4", tmp2, "1", "0", "2")
    # l_("5") = 3maj(3maj(c, mdef, b), mdef, a)
    hdc_ternary_majority_const_vr(tmp2, "5", tmp2, "4", "0", "3")

    hdc_cp(tmp2, "6", i, s)
    hdc_cp(tmp2, "7", h, s)
    hdc_cp(tmp2, "8", g, s)
    # l___("9") = 3maj(i, l_, h)
    hdc_ternary_majority_const_vr(tmp2, "9", tmp2, "6", "5", "7")
    # l("A") = 3maj(3maj(i, l_, h), l_, g)
    hdc_ternary_majority_const_vr(tmp2, "A", tmp2, "9", "5", "8")

    # We're done with t2mdef, i, h, g, l_, l__, and l___, so
    # 0, 4, 5, 6, 7, 8, 9 are dead (free)
    hdc_cp(tmp2, "0", tmp1, t1mghi)
    # r__("4") = 3maj(c, mghi, b)
    hdc_ternary_majority_const_vr(tmp2, "4", tmp2, "1", "0", "2")
    # r_("5") = 3maj(3maj(c, mghi, b), mghi, a)
    hdc_ternary_majority_const_vr(tmp2, "5", tmp2, "4", "0", "3")

    hdc_cp(tmp2, "6", f, s)
    hdc_cp(tmp2, "7", e, s)
    hdc_cp(tmp2, "8", d, s)
    # r___("9") = 3maj(f, r_, e)
    hdc_ternary_majority_const_vr(tmp2, "9", tmp2, "6", "5", "7")
    # r("B") = 3maj(3maj(f, r_, e), r_, d)
    hdc_ternary_majority_const_vr(tmp2, "B", tmp2, "9", "5", "8")

    t1m : Section = Belex.Section(3)
    hdc_ternary_majority_const_vr(
        tmp1, t1m, tmp1, t1mabc, t1mdef, t1mghi)

    t1l: Section = Belex.Section(4)
    hdc_cp(tmp1, t1l, tmp2, "A")

    t1r: Section = Belex.Section(5)
    hdc_cp(tmp1, t1r, tmp2, "B")

    hdc_ternary_majority_const_vr(
        dv, s, tmp1, t1l, t1m, t1r)

    pass


@belex_apl
def hdc_9_maj_via_3_const_vr(
        Belex,
        dv   : VR, ds   : Section,
        tmp1 : VR, tmp2 : VR) -> None:
    """Call precalc in the caller, then call this."""

    t1mabc : Section = Belex.Section(0)
    t1mdef : Section = Belex.Section(1)
    t1mghi : Section = Belex.Section(2)

    t1m : Section = Belex.Section(3)
    hdc_ternary_majority_const_vr(
        tmp1, t1m, tmp1, t1mabc, t1mdef, t1mghi)

    t1l : Section = Belex.Section(4)
    hdc_cp(tmp1, t1l, tmp2, "A")

    t1r : Section = Belex.Section(5)
    hdc_cp(tmp1, t1r, tmp2, "B")

    hdc_ternary_majority_const_vr(
        dv, ds, tmp1, t1l, t1m, t1r)


@belex_apl
def hdc_9_maj_via_3_const_vr_precalc(
        Belex,
        sv : VR,
        a : Section, b : Section, c : Section,
        d : Section, e : Section, f : Section,
        g : Section, h : Section, i : Section,
        tmp1: VR, tmp2: VR) -> None:
    """Leaves info in known sections, 0, 1, 2."""

    # import pudb; pudb.set_trace();

    t1mabc : Section = Belex.Section(0)
    hdc_ternary_majority_const_vr(tmp1, t1mabc, sv, a, b, c)

    t1mdef: Section = Belex.Section(1)
    hdc_ternary_majority_const_vr(tmp1, t1mdef, sv, d, e, f)

    t1mghi: Section = Belex.Section(2)
    hdc_ternary_majority_const_vr(tmp1, t1mghi, sv, g, h, i)

    # String section numbers as hex digits are in tmp2.
    # Notice that the roles of VR and Section are reversed compared
    # to the parallel-looking code in "hdc_9_maj_via_3_const_section."
    hdc_cp(tmp2, "0", tmp1, t1mdef)
    hdc_cp(tmp2, "1", sv, c)
    hdc_cp(tmp2, "2", sv, b)
    hdc_cp(tmp2, "3", sv, a)
    # l__("4") = 3maj(c, mdef, b)
    hdc_ternary_majority_const_vr(tmp2, "4", tmp2, "1", "0", "2")
    # l_("5") = 3maj(3maj(c, mdef, b), mdef, a)
    hdc_ternary_majority_const_vr(tmp2, "5", tmp2, "4", "0", "3")

    hdc_cp(tmp2, "6", sv, i)
    hdc_cp(tmp2, "7", sv, h)
    hdc_cp(tmp2, "8", sv, g)
    # l___("9") = 3maj(i, l_, h)
    hdc_ternary_majority_const_vr(tmp2, "9", tmp2, "6", "5", "7")
    # l("A") = 3maj(3maj(i, l_, h), l_, g)
    hdc_ternary_majority_const_vr(tmp2, "A", tmp2, "9", "5", "8")

    # We're done with t2mdef, i, h, g, l_, l__, and l___, so
    # 0, 4, 5, 6, 7, 8, 9 are dead (free)
    hdc_cp(tmp2, "0", tmp1, t1mghi)
    # r__("4") = 3maj(c, mghi, b)
    hdc_ternary_majority_const_vr(tmp2, "4", tmp2, "1", "0", "2")
    # r_("5") = 3maj(3maj(c, mghi, b), mghi, a)
    hdc_ternary_majority_const_vr(tmp2, "5", tmp2, "4", "0", "3")

    hdc_cp(tmp2, "6", sv, f)
    hdc_cp(tmp2, "7", sv, e)
    hdc_cp(tmp2, "8", sv, d)
    # r___("9") = 3maj(f, r_, e)
    hdc_ternary_majority_const_vr(tmp2, "9", tmp2, "6", "5", "7")
    # r("B") = 3maj(3maj(f, r_, e), r_, d)
    hdc_ternary_majority_const_vr(tmp2, "B", tmp2, "9", "5", "8")

    # About to exhaust SM registers. Make a new fragment.

    # t1m : Section = Belex.Section(3)
    # hdc_ternary_majority_const_vr(
    #     tmp1, t1m, tmp1, t1mabc, t1mdef, t1mghi)

    # t1l: Section = Belex.Section(4)
    # hdc_cp(tmp1, t1l, tmp2, "A")
    #
    # t1r: Section = Belex.Section(5)
    # hdc_cp(tmp1, t1r, tmp2, "B")

    return


#   _ ___     __  __       _         _ _
#  / | __|___|  \/  |__ _ (_)___ _ _(_) |_ _  _
#  | |__ \___| |\/| / _` || / _ \ '_| |  _| || |
#  |_|___/   |_|  |_\__,_|/ \___/_| |_|\__|\_, |
#                       |__/               |__/


@belex_apl
def hdc_15_maj_const_vr(
        Belex,
        dest_vr: VR, dest_section: Section,
        scratch_vr: VR,
        av_vr: VR
) -> None:
    dest_vr[::] <= RSP16()
    scratch_vr[::] <= RSP16()

    RL["0"] <= av_vr()
    GL["0"] <= RL()
    dest_vr["0"] <= GL()

    glass_kwargs = {'plats': 30, 'sections': 16,
                    'fmt': 'hex', 'order': 'lsb'}

    for sec in ["1", "2", "3", "4", "5", "6", "7",
                "8", "9", "A", "B", "C", "D", "E"]:
        RL[sec] <= av_vr()
        GL[sec] <= RL()
        scratch_vr["0"] <= GL()
        add_u16(dest_vr, dest_vr, scratch_vr)

    # If bit 3 is ON, the count is 8 or greater, i.e., a majority
    cpy_imm_16(scratch_vr, 0x0003)

    Belex.glass(scratch_vr, **glass_kwargs)
    Belex.glass(dest_vr, **glass_kwargs)

    RL[::] <= dest_vr()
    RL[::] &= scratch_vr()
    GL[0x0003] <= RL()

    Belex.glass(RL, **glass_kwargs)

    dest_vr[dest_section] <= GL()
    pass


@belex_apl
def hdc_15_maj_const_vr_first_eight_sections(
        Belex,
        dest_vr: VR, dest_section: Section,
        scratch_vr: VR,
        av_vr: VR
) -> None:
    """Break up 15-maj so we don't run out of section reg's."""
    dest_vr[::] <= RSP16()
    scratch_vr[::] <= RSP16()

    RL["0"] <= av_vr()
    GL["0"] <= RL()
    dest_vr["0"] <= GL()

    for sec in ["1", "2", "3", "4", "5", "6", "7"]:
        RL[sec] <= av_vr()
        GL[sec] <= RL()
        scratch_vr["0"] <= GL()
        add_u16(dest_vr, dest_vr, scratch_vr)

    pass


@belex_apl
def hdc_15_maj_const_vr_last_seven_sections(
        Belex,
        dest_vr: VR, dest_section: Section,
        scratch_vr: VR,
        av_vr: VR
) -> None:
    """Break up 15-maj so we don't run out of section reg's."""
    for sec in ["8", "9", "A", "B", "C", "D", "E"]:
        RL[sec] <= av_vr()
        GL[sec] <= RL()
        scratch_vr["0"] <= GL()
        add_u16(dest_vr, dest_vr, scratch_vr)

    # If bit 3 is ON, the count is 8 or greater, i.e., a majority
    cpy_imm_16(scratch_vr, 0x0003)

    RL[::] <= dest_vr()
    RL[::] &= scratch_vr()
    GL[0x0003] <= RL()

    dest_vr[dest_section] <= GL()
    pass


@belex_apl
def hdc_15_maj_const_vr_eight_explicit_sections(
        Belex,
        dest_vr: VR, dest_section: Section,
        scratch_vr: VR,
        av_vr: VR,
        a: Section, b: Section, c: Section, d: Section,
        e: Section, f: Section, g: Section, h: Section,
) -> None:
    """Break up 15-maj so we don't run out of section reg's."""
    dest_vr[::] <= RSP16()
    scratch_vr[::] <= RSP16()

    RL[a] <= av_vr()
    GL[a] <= RL()
    dest_vr["0"] <= GL()

    for sec in [b, c, d, e, f, g, h]:
        RL[sec] <= av_vr()
        GL[sec] <= RL()
        scratch_vr["0"] <= GL()
        add_u16(dest_vr, dest_vr, scratch_vr)

    pass


@belex_apl
def hdc_15_maj_const_vr_seven_explicit_sections(
        Belex,
        dest_vr: VR, dest_section: Section,
        scratch_vr: VR,
        av_vr: VR,
        i: Section, j: Section, k: Section, l: Section,
        m: Section, n: Section, o: Section
) -> None:
    """Break up 15-maj so we don't run out of section reg's."""
    for sec in [i, j, k, l, m, n, o]:
        RL[sec] <= av_vr()
        GL[sec] <= RL()
        scratch_vr["0"] <= GL()
        add_u16(dest_vr, dest_vr, scratch_vr)

    # If bit 3 is ON, the count is 8 or greater, i.e., a majority
    cpy_imm_16(scratch_vr, 0x0003)

    RL[::] <= dest_vr()
    RL[::] &= scratch_vr()
    GL[0x0003] <= RL()

    dest_vr[dest_section] <= GL()
    pass


#   ___                _   __  __      _ _   _      _           _      _    _
#  | __|  _ ___ ___ __| | |  \/  |_  _| | |_(_)_ __| |_  _ ___ /_\  __| |__| |
#  | _| || (_-</ -_) _` | | |\/| | || | |  _| | '_ \ | || |___/ _ \/ _` / _` |
#  |_| \_,_/__/\___\__,_| |_|  |_|\_,_|_|\__|_| .__/_|\_, |  /_/ \_\__,_\__,_|
#                                             |_|     |__/


@belex_apl
def hdc_ternary_fma_const_vr (
        Belex,
        dst : VR, d : Section,  # also scratch in src
        src : VR, a : Section,
                  b : Section,
                  c : Section,
                  x : Section) -> None:
    """Compute x * (a + b + c), where each bhv
    inhabits a section of src VR. Deposit result is
    dst[d], where d does not equal any of the other
    sections. src[d] is overwritten as scratch
    space. Exploit function inlining in Belex.
    """
    hdc_ternary_majority_const_vr(dst, d, src, a, b, c)
    hdc_cp(src, d, dst, d)
    hdc_mul_const_vr(dst, d, src, x, d)


@belex_apl
def hdc_ternary_fma_const_section(
        Belex,
        dst : VR,
        a   : VR,
        b   : VR,
        c   : VR,
        x   : VR,
        sec : Section) -> None:
    """Compute x * (a + b + c), where each bhv
    inhabits a section of src VR. Deposit result is
    dst[d], where d does not equal any of the other
    sections. src[d] is overwritten as scratch
    space. Exploit function inlining in Belex.
    """
    hdc_ternary_majority_const_section(dst, a, b, c, sec)
    hdc_mul_const_section(dst, dst, x, sec)


#   _  _                  _                 _
#  | || |__ _ _ __  _ __ (_)_ _  __ _   ___| |_ __
#  | __ / _` | '  \| '  \| | ' \/ _` | / -_)  _/ _|_
#  |_||_\__,_|_|_|_|_|_|_|_|_||_\__, | \___|\__\__(_)
#                               |___/


@belex_apl
def hdc_hamming_henry_warren(
        Belex,
        dst : VR,
        # inputs
        x0 : VR, s0 : Section,
        # combs; stored in sections 0 through A
        c: VR,
) -> None:
    """
    UNDONE:

    Method of combs without nocks. With a 32-bit word: notice
    the '+'
               c0...cA precomputed
    x1 = (x0 & 0x55555555) + ((x0 >>    1) & 0x55555555);
    x2 = (x1 & 0x33333333) + ((x1 >>    2) & 0x33333333);
    x3 = (x2 & 0x0F0F0F0F) + ((x2 >>    4) & 0x0F0F0F0F);
    x4 = (x3 & 0x00FF00FF) + ((x3 >>    8) & 0x00FF00FF);
    x5 = (x4 & 0x0000FFFF) + ((x4 >>   16) & 0x0000FFFF);

    We work with a "word" of 2048 bits.
    x6 = (x5 & 0x00....FF) + ((x5 >>   32) & 0x00....FF);
    x7 = (x6 & 0x00....FF) + ((x6 >>   64) & 0x00....FF);
    x8 = (x7 & 0x00....FF) + ((x7 >>  128) & 0x00....FF);
    x9 = (x8 & 0x00....FF) + ((x8 >>  256) & 0x00....FF);
    xA = (x9 & 0x00....FF) + ((x9 >>  512) & 0x00....FF);
    xB = (xA & 0x00....FF) + ((xA >> 1024) & 0x00....FF);

    Return the accumulator, containing xB, to the device
    controller, which will do final sums of the 2048-bit
    word counts into 8192-bit superword counts and into
    32K-bit hyperword counts.
    """

    x = Belex.VR(0)  # accumulator
    y = Belex.VR(0)  # shifted

    # clear
    x[:] <= RSP16()
    y[:] <= RSP16()
    dst[:] <= RSP16()

    RL["0"] <= c()
    GL["0"] <= RL()        # 0x55555555...
    RL[s0] <= x0() & GL()  # (x0 & 0x55555555...)
    GL[s0] <= RL()
    x["0"] <= GL()

    RL[s0] <= x0()
    y[s0] <= WRL()
    RL[s0] <= y() & GL()  # ((x0 >> 1) & 0x55555555...)
    GL[s0] <= RL()
    y["0"] <= GL()

    add_u16(dst, x, y)

    RL["1"] <= c()
    GL["1"] <= RL()        # 0x33333333...
    RL[:] <= dst() & GL()  # (x1 & 0x33333333...)
    GL[s0] <= RL()
    x["0"] <= GL()

    RL[s0] <= x0()
    y[s0] <= WRL()
    RL[s0] <= y() & GL()  # ((x1 >> 1) & 0x33333333...)
    GL[s0] <= RL()
    y["0"] <= GL()


@belex_apl
def hdc_hamming_prepare_input_for_writeback(
        Belex,
        input_vr: VR, input_sect: Section
) -> None:
    # Clear junk in the input_vr for future writeback of counts.
    input_vr[~input_sect] <= RSP16()


@belex_apl
def hdc_combs_and_nocks_8x(
        Belex,
        c0: VR, c1: VR, c2: VR, c3: VR,
        c4: VR, c5: VR, c6: VR, c7: VR,
        # c8: VR, c9: VR, cA: VR, cB: VR,
        # cC: VR, cD: VR, cE: VR, cF: VR,
        input_vr: VR,  input_sect: Section,
        resonators : Section,  # section 0 of multiple VRs
        nocks      : Section,  # section 1 of multiple VRs
        comb       : VR,
        # bit counts: 4 bits, little-endian, cross-section
        bc0 : Section, bc1 : Section,
        bc2 : Section, bc3 : Section,) -> Section:

    # Let r = a resonator, that is, a section full of repeating
    # 4-bit patterns, e.g., 0101 0101 0101 ... = 5 5 5 ... = c5 in
    # big-endian, left-to-right. Let i be an input vector. Think
    # of it as divided into big-endian nibbles, e.g., 1101 0101
    # 1001 ... = D 5 9 ... = input_vr[input_sect]. For each nibble
    # k, compute e[k] = r[k] XOR i[k], which is 0 iff r[k] ==
    # i[k]. Detect 0 by OR'ing up the bits:
    # z[k] = e[k] | (e[k] << 1) | (e[k] << 2) | (e[k] << 3).
    # Produce counts by looking them up little-endian in the bc
    # sections of the combs VRs.

    glass_kwargs = {'orientation': 'section-wise', 'fmt':'hex',
                    'plats': 120, 'sections': 16, 'order': 'msb'}

    # 'input_glass' in the test:
    Belex.glass(input_vr, **glass_kwargs)

    def resonate_comb_and_nock(resonator: VR) -> None:
        RL[input_sect] <= input_vr()
        GL[input_sect] <= RL()

        RL[resonators]   <= resonator() ^ GL()
        GL[resonators]   <= RL()

        RL[nocks]        <= GL()
        resonator[nocks] <= RL()
        resonator[nocks] |= ERL()
        RL[nocks]        <= resonator() & comb()

        RL[nocks]        <= resonator()
        resonator[nocks] |= ERL()
        RL[nocks]        <= resonator() & comb()

        RL[nocks]        <= resonator()
        resonator[nocks] |= ERL()
        RL[nocks]        <= resonator() & comb()

        # Spread the nock bit back out Eest, overwriting
        # garbage left there by the ERLs above.
        RL[nocks]        |= WRL()
        RL[nocks]        |= WRL()
        RL[nocks]        |= WRL()

        resonator[nocks] <= RL()
        GL[nocks]        <= RL()
        # workaround Issue 79 -- https://bitbucket.org/gsitech/belex/issues/79/index-calc-works-in-pytest-but-not-in
        bit_count_sections = "2345"
        RL[bit_count_sections] <= resonator() & ~GL()
        resonator[bit_count_sections] <= RL()
        input_vr[bit_count_sections] |= RL()

    # Test on c5, trust the rest later:

    resonate_comb_and_nock(resonator=c5)

    # 'c5_glass' in the test:
    # Should see 0 exactly where input == 5
    Belex.glass(c5, **glass_kwargs)

    resonate_comb_and_nock(resonator=c0)
    resonate_comb_and_nock(resonator=c1)
    # Debug c1:
    Belex.glass(c1, **glass_kwargs)

    resonate_comb_and_nock(resonator=c2)
    resonate_comb_and_nock(resonator=c3)
    resonate_comb_and_nock(resonator=c4)
    # c5 was already done.
    resonate_comb_and_nock(resonator=c6)
    resonate_comb_and_nock(resonator=c7)

    # Do the input again to see the bit-counts written back
    Belex.glass(input_vr, **glass_kwargs)

    return nocks


@belex_apl
def hdc_hamming_sum(Belex,
                    input_vr: VR, output_vr: VR, temp_vr: VR,
                    # sections implicit due to Issue #79
                    ) -> None:
    sum_sections = "23456789ABCDEF"
    def summup(into_: VR, from_: VR, sumd_: VR,
               n_plus_1: int):
        # Clear temp and output:
        temp_vr[:] <= RSP16()
        RL[sum_sections] <= from_()
        temp_vr[sum_sections] <= ERL()
        for i in range(n_plus_1 - 1):
            RL[sum_sections] <= temp_vr()
            temp_vr[sum_sections] <= ERL()
        add_u16(into_, sumd_, temp_vr)

    output_vr[:] <= RSP16()

    summup(output_vr, input_vr, input_vr, 4)
    summup(output_vr, input_vr, output_vr, 8)
    summup(output_vr, input_vr, output_vr, 12)

    summup(output_vr, output_vr, output_vr, 16)
    summup(output_vr, output_vr, output_vr, 32)
    summup(output_vr, output_vr, output_vr, 64)
    summup(output_vr, output_vr, output_vr, 128)
    summup(output_vr, output_vr, output_vr, 256)
    summup(output_vr, output_vr, output_vr, 512)
    pass


@belex_apl
def hdc_hamming_sum_tuning(Belex,
                    input_vr: VR, output_vr: VR, temp_vr: VR,
                    # sections implicit due to Issue #79
                    ) -> None:
    """Experiment with reducing clocks in the MMB and promoting
    more final sums to the ARC."""
    sum_sections = "23456789ABCDEF"
    def summup(into_: VR, from_: VR, sumd_: VR,
               n_plus_1: int):
        # Clear temp and output:
        temp_vr[:] <= RSP16()
        RL[sum_sections] <= from_()
        temp_vr[sum_sections] <= ERL()
        for i in range(n_plus_1 - 1):
            RL[sum_sections] <= temp_vr()
            temp_vr[sum_sections] <= ERL()
        add_u16(into_, sumd_, temp_vr)

    output_vr[:] <= RSP16()

    summup(output_vr, input_vr, input_vr, 4)
    summup(output_vr, input_vr, output_vr, 8)
    summup(output_vr, input_vr, output_vr, 12)

    summup(output_vr, output_vr, output_vr, 16)
    summup(output_vr, output_vr, output_vr, 32)
    summup(output_vr, output_vr, output_vr, 64)
    summup(output_vr, output_vr, output_vr, 128)
    # summup(output_vr, output_vr, output_vr, 256)
    # summup(output_vr, output_vr, output_vr, 512)
    pass


#   ___                   _        _   _
#  | _ \___ _ _ _ __ _  _| |_ __ _| |_(_)___ _ _  ___
#  |  _/ -_) '_| '  \ || |  _/ _` |  _| / _ \ ' \(_-<
#  |_| \___|_| |_|_|_\_,_|\__\__,_|\__|_\___/_||_/__/


# A BHV, aka _vector_, is always of length 8192 bits. Think of a permutation
# vector, P, as containing, horizontally, 1024 8‑bit bytes, P[i], for all
# i ∈ {0, 1, …, 1023}, such that each integer between 0 and 1023, inclusive
# both ends, occurs exactly once in P. _Applying_ a permutation P to a vector
# A produces a vector, B, whose i‑th element contains A[P[i]]:
#
#     B = P(A) <=> B[i] = A[P[i]] ∀ i ∈ {0, 1, …, 1023}                    (1)
#
# Define the inverse permutation, ~P, of P, so that
#
#     A = ~P(B) <=> A[i] = B[~P[i]] ∀ i ∈ {0, 1, …, 1023}                  (2)
#
# Let P[i] = j.
#
#     B[i] = A[j] = B[~P[j]]                                               (3)
#
# because Equation 2 holds for any i, for j in particular. Equation 3 implies
#
#     ~P[j] = ~P[P[i]] = i                                                 (4)
#
# which is a construction for ~P.
#
# To illustrate, consider vectors of 3-bytes. There are six permutations:
#
#     P0 = [0, 1, 2]  # The identity permutation
#     P1 = [0, 2, 1]
#     P2 = [1, 0, 2]
#     P3 = [1, 2, 0]
#     P4 = [2, 0, 1]
#     P5 = [2, 1, 0]
#
# Choose P = P2 for discussion.
#
#       A      P = P2     B, Equation 1
#     +------+----------+-----------------------+
#     | A[0] | P[0] = 1 | B[0] = A[P[0]] = A[1] |
#     | A[1] | P[1] = 0 | B[1] = A[P[1]] = A[0] |
#     | A[2] | P[2] = 2 | B[2] = A[P[2]] = A[2] |
#     +------+----------+-----------------------+
#
#       ~P = ~P2, Equation 4   Equation 2
#     +----------------------+------------------------+
#     | ~P[0] = ~P[P[1]] = 1 | B[~P[0]] = B[1] = A[0] |
#     | ~P[1] = ~P[P[0]] = 0 | B[~P[1]] = B[0] = A[1] |
#     | ~P[2] = ~P[P[2]] = 2 | B[~P[2]] = B[2] = A[2] |
#     +----------------------+------------------------+
#
# In this case, P2 == ~P2, an impertinent fact due to P's having just one
# 2-cycle, a _swap_. Now choose P = P4:
#
#       A      P = P4     B, Equation 1
#     +------+----------+-----------------------+
#     | A[0] | P[0] = 2 | B[0] = A[P[0]] = A[2] |
#     | A[1] | P[1] = 0 | B[1] = A[P[1]] = A[0] |
#     | A[2] | P[2] = 1 | B[2] = A[P[2]] = A[1] |
#     +------+----------+-----------------------+
#
#       ~P = ~P4, Equation 4   Equation 2
#     +----------------------+------------------------+
#     | ~P[0] = ~P[P[1]] = 1 | B[~P[0]] = B[1] = A[0] |
#     | ~P[1] = ~P[P[2]] = 2 | B[~P[1]] = B[2] = A[1] |
#     | ~P[2] = ~P[P[0]] = 0 | B[~P[2]] = B[0] = A[2] |
#     +----------------------+------------------------+
#
# In this case, P4 = [2, 0, 1] and ~P4 = [1, 2, 0], with no trivial
# relationship.
#
# A similar analysis hold if permutations are taken not byte-wise, but
# bit-wise, with 8192 indices instead of 1024. The difference is that the
# permutation indices cannot be stored _horizontally_ in an 8‑Kibit vector.
# Instead, in the APU, we may store the indices _vertically_, across 13
# sections.

@belex_apl
def hdc_apply_permutation(
        Belex,
        dest_vr : VR, dest_sc : Section,
        perm_vr : VR, perm_sc : Section,
        vict_vr : VR, vict_se : Section,
) -> None:
    pass
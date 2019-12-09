// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2019, Sara Hetzel <hetzel @ molgen.mpg.de>
// Copyright (c) 2016-2019, Knut Reinert and Freie Universität Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is Free Software: you can redistribute it and/or modify it
// under the terms found in the LICENSE[.md|.rst] file distributed
// together with this file.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// ==========================================================================
// Scoring schemes for bisulfite converted data on dna5 alphabet
// ==========================================================================

#pragma once

#include <range/v3/algorithm/copy.hpp>

#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/algorithm>

namespace seqan3
{

template <Arithmetic score_type = int8_t>
class bisulfite_scoring_scheme : public scoring_scheme_base<bisulfite_scoring_scheme<score_type>, dna5, score_type>
{
private:
    using base_t = scoring_scheme_base<bisulfite_scoring_scheme<score_type>, dna5, score_type>;
    using base_t::matrix;
public:
    using base_t::base_t;
    using typename base_t::matrix_type;
    using matrix_size_type = std::remove_const_t<decltype(alphabet_size<dna5>)>;
    static constexpr matrix_size_type matrix_size = alphabet_size<dna5>;

    constexpr bisulfite_scoring_scheme() noexcept {}

    template <Arithmetic score_arg_t>
    constexpr bisulfite_scoring_scheme(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    {
        set_bisulfite_scheme(ms, mms);
    }

    constexpr bisulfite_scoring_scheme(matrix_type const & _matrix) noexcept {}

    template <Arithmetic score_arg_t>
    constexpr void set_bisulfite_scheme(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    {
        std::conditional_t<std::Integral<score_type>, int64_t, double> i_ms = static_cast<score_arg_t>(ms);
        std::conditional_t<std::Integral<score_type>, int64_t, double> i_mms = static_cast<score_arg_t>(mms);
        if ((i_ms  < std::numeric_limits<score_type>::lowest() || i_ms  > std::numeric_limits<score_type>::max()) ||
            (i_mms < std::numeric_limits<score_type>::lowest() || i_mms > std::numeric_limits<score_type>::max()))
        {
            throw std::invalid_argument{"You passed a score value to set_bisulfite_scheme that is out of range of the "
                                        "scoring scheme's underlying type. Define your scoring scheme with a larger "
                                        "template parameter or down-cast you score value beforehand to prevent "
                                        "this exception."};
        }

        // Assume query is horizontal sequence
        for (matrix_size_type i = 0; i < matrix_size; ++i)
            for (matrix_size_type j = 0; j < matrix_size; ++j)
                matrix[i][j] = (i == j || (i == 4 && j == 1)) ? static_cast<score_type>(i_ms) : static_cast<score_type>(i_mms);
    }
};

bisulfite_scoring_scheme() -> bisulfite_scoring_scheme<int8_t>;

template <Arithmetic score_arg_type>
bisulfite_scoring_scheme(match_score<score_arg_type>,
                         mismatch_score<score_arg_type>) -> bisulfite_scoring_scheme<int8_t>;

template <Arithmetic score_arg_type>
bisulfite_scoring_scheme(std::array<std::array<score_arg_type, 5>, 5>) -> bisulfite_scoring_scheme<score_arg_type>;

} // namespace seqan3

namespace seqan
{

struct Bisulfite{};

template <typename TValue>
class Score<TValue, Bisulfite> {
public:
    // The score for a match.
    TValue data_match;

    // The score for a mismatch.
    TValue data_mismatch;

    // The gap extension score.
    TValue data_gap_extend;

    // The gap open score.
    TValue data_gap_open;

    Score()
        : data_match(0), data_mismatch(-1), data_gap_extend(-1),
          data_gap_open(-1) {
    }

    Score(TValue _match, TValue _mismatch, TValue _gap)
        : data_match(_match), data_mismatch(_mismatch),
          data_gap_extend(_gap), data_gap_open(_gap) {
    }

    Score(TValue _match, TValue _mismatch, TValue _gap_extend, TValue _gap_open)
        : data_match(_match), data_mismatch(_mismatch),
          data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
    }
};

typedef Score<int, Bisulfite> BisulfiteScore;

inline int
scoreMatch(BlastScoringScheme<Score<int, Bisulfite>> & scheme)
{
    return scoreMatch(scheme._internalScheme);
}

inline int
scoreMismatch(BlastScoringScheme<Score<int, Bisulfite>> & scheme)
{
    return scoreMismatch(scheme._internalScheme);
}

inline bool
_selectSet(BlastScoringScheme<Score<int, Bisulfite>> & scheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;

    for (typename TKAValues::TSize i = 0; i < TKAValues::nParamSets; ++i)
    {
        if ((TKAValues::VALUE[i][0] ==  scoreMatch(scheme))  &&
            (TKAValues::VALUE[i][1] == -scoreMismatch(scheme)) &&
            (TKAValues::VALUE[i][2] == -scoreGapOpenBlast(scheme)) &&
            (TKAValues::VALUE[i][3] == -scoreGapExtend(scheme)))
        {
            scheme.parameterIndex = i;
            return true;
        }
    }
    // no suitable adapter
    scheme.parameterIndex = std::numeric_limits<typename TKAValues::TSize>::max();
    return false;
}

inline bool
setScoreMatch(BlastScoringScheme<Score<int, Bisulfite>> & scheme, int const val)
{
    setScoreMatch(scheme._internalScheme, val);
    return _selectSet(scheme);
}

inline bool
setScoreMismatch(BlastScoringScheme<Score<int, Bisulfite>> & scheme, int const val)
{
    setScoreMismatch(scheme._internalScheme, val);
    return _selectSet(scheme);
}

inline double
getLambda(BlastScoringScheme<Score<int, Bisulfite>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][4];
}

inline double
getKappa(BlastScoringScheme<Score<int, Bisulfite>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][5];
}

inline double
getH(BlastScoringScheme<Score<int, Bisulfite>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Bisulfite>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][6];
}

inline double
getAlpha(BlastScoringScheme<Score<int, Bisulfite>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][7];
}

inline double
getBeta(BlastScoringScheme<Score<int, Bisulfite>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][8];
}

template <typename TValue, typename TSpec, typename TSeqHVal, typename TSeqVVal>
    requires std::Same<TSpec, Bisulfite>
inline TValue
score(Score<TValue, TSpec> const & me, TSeqHVal valH, TSeqVVal valV) {
    if (valH == valV || (valH == 'T' && valV == 'C'))
        return scoreMatch(me);
    else
        return scoreMismatch(me);
}

} // namespace seqan

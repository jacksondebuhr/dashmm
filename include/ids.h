#ifndef __DASHMM_BUILTIN_IDENTIFIERS_H__
#define __DASHMM_BUILTIN_IDENTIFIERS_H__


/// \file include/builtin_ids.h
/// \brief Declarations of built-in expansion and method identifiers


namespace dashmm {


/// The first built-in method type identifier
constexpr int kFirstMethodType = 0;

/// The last built-in method type identifier
constexpr int kLastMethodType = 999;

/// The first user method type identifier
constexpr int kFirstUserMethodType = 1000;

/// The last user method type identifier
constexpr int kLastUserMethodType = 1999;

/// The identifier for the BH Method
constexpr int kMethodBH = kFirstMethodType;
// etc...


/// The lowest allowed expansion identifier for internal use
constexpr int kFirstExpansionType = 2000;

/// The highest allowed expansion identifier for internal use
constexpr int kLastExpansionType = 2999;

/// The first allowed expansion identifier for user defined expansions
constexpr int kFirstUserExpansionType = 3000;

/// the highest allowed expansion identifier for user defined expansions
constexpr int kLastUserExpansionType = 3999;

/// The identifier for the Laplace Center of Mass Expansion
constexpr int kExpansionLaplaceCOM = kFirstExpansionType;
//constexpr int kExpansionEtc.. = kExpansionLaplaceCOM + 1;



} // namespace dashmm


#endif // __DASHMM_BUILTIN_IDENTIFIERS_H__

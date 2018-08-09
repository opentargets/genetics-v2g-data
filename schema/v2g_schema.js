type Query {
    v2g(variantId: String!, evidenceFilter: [String!], returnEvidence: Bool!): VariantAggregate
    // i2g(indexVariantId: String!, studyId: String!, evidenceFilter: [String!], returnEvidence: Bool!): IndexAggregate
}

/*
// Per indexVariant assignment
type IndexAggregate {
    indexVariantId: String!
    studyId: String!
    iaRank: [GeneRankScore!]
    iaEvidence: [VariantAggregate!]
}
*/

// Per variant assignment
type VariantAggregate {
    variantId: String!

    vaRank: [GeneRankScore!]
    vaEvidence: [SourceAggregate!]
}

// Per source assignment
type SourceAggregate {
    variantId: String!
    sourceId: String!

    saRank: [GeneRankScore!]
    saEvidence: [TissueEvidence!]
}

// Per tissue evidence
type TissueEvidence {
    variantId: String!
    sourceId: String!
    tissueId: String!

    rawEvidence: [GeneRawScore!]
}

// Object containing score produced by ranking algorithm
type GeneRankScore {
    geneName: String!
    geneId: String!

    score: Float! // Number 0-1 used for ranking but will be displayed as [High|Medium|Low] evidence
}

// Object containing raw experiment score / quantile
type GeneRawScore {
    variantId: String!
    sourceId: String!
    tissueId: String!
    geneName: String!
    geneId: String!

    quantile: Float!
    beta: Float // Only for QTL evidence (eQTL, pQTL)
    pval: Float // Only for QTL evidence (eQTL, pQTL)
    intervalScore: Float // Only for Interval evidence (PCHiC, Fantom5, etc.)
}

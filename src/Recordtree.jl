abstract type AbstractNode end

mutable struct Directive <: AbstractNode
    record::GFF3.Record
end

AbstractTrees.children(::Directive) = ()

mutable struct Comment <: AbstractNode
    record::GFF3.Record
end

AbstractTrees.children(::Comment) = ()

mutable struct Feature <: AbstractNode
    record::GFF3.Record
    seqid::Union{Missing, String}
    featuretype::Union{Missing, String}
    source::Union{Missing, String}
    seqstart::Union{Missing, Int}
    seqend::Union{Missing, Int}
    score::Union{Missing, Float64}
    strand::Union{Missing, GFF3.GenomicFeatures.Strand}
    phase::Union{Missing, Int}
    attributes::Vector{Pair{String,Vector{String}}}
    children::Vector{Feature}
end

AbstractTrees.children(f::Feature) = f.children

mutable struct Chromosome <: AbstractNode
    record::String
    children::Vector{Feature}
end
AbstractTrees.children(c::Chromosome) = c.children


AbstractTrees.nodevalue(n::AbstractNode) = n.record

AbstractTrees.printnode(io::IO, d::Directive) = print(io, d.record)
AbstractTrees.printnode(io::IO, c::Comment) = print(io, c.record)
AbstractTrees.printnode(io::IO, f::Feature) = print(io, f.record)
AbstractTrees.printnode(io::IO, c::Chromosome) = print(io, c.record)

function Chromosome(c::String)
    children = Vector{Feature}[]
    return Chromosome(c, children)
end

function Feature(record::GFF3.Record)
    seqid = GFF3.hasseqid(record) ? GFF3.seqid(record) : missing
    featuretype = GFF3.hasfeaturetype(record) ? GFF3.featuretype(record) : missing
    source = GFF3.hassource(record) ? GFF3.source(record) : missing
    seqstart = GFF3.hasseqstart(record) ? GFF3.seqstart(record) : missing
    seqend = GFF3.hasseqend(record) ? GFF3.seqend(record) : missing
    score = GFF3.hasscore(record) ? GFF3.score(record) : missing
    strand = GFF3.hasstrand(record) ? GFF3.strand(record) : missing
    phase = GFF3.hasphase(record) ? GFF3.phase(record) : missing
    attributes = GFF3.attributes(record)
    children = Vector{Feature}[]
    return Feature(record, seqid, featuretype, source, seqstart, seqend, score, strand,  phase, attributes, children)
end

function update!(f::Feature)
    attributes_str = join(map(attr -> first(attr) * "=" * join(map(x -> URIParser.escape(x), last(attr)), ","), f.attributes), ";")
    seqid = encode(f.seqid)
    featuretype = encode(f.featuretype)
    source = encode(f.source)
    seqstart = encode(f.seqstart)
    seqend = encode(f.seqend)
    score = encode(f.score)
    strand = encode(f.strand)
    phase = encode(f.phase)
    attributes = encode(attributes_str)

    content = join(
        [seqid,
        source,
        featuretype,
        seqstart,
        seqend,
        score,
        strand,
        phase,
        attributes
        ],
        "\t"
    )
    @show content
    new_record = GFF3.Record(content)
    f.record = new_record
end

function add_children!(c::Chromosome, children::Vector{Feature})
    append!(c.children,children)
end

function add_children!(f::Feature, children::Vector{Feature})
    append!(f.children,children)
end

encode(m::Missing) = "."
encode(s::String) = isempty(s) ? "." : s
encode(n::Number) = string(n) 
function encode(g::GFF3.GenomicFeatures.Strand)
    if g == GFF3.GenomicFeatures.STRAND_NA
        return "?"
    elseif g == GFF3.GenomicFeatures.STRAND_POS
        return "+"
    elseif g == GFF3.GenomicFeatures.STRAND_NEG
        return "-"
    elseif g ==GFF3.GenomicFeatures.STRAND_BOTH
        return "."
    end
end
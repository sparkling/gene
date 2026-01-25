<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  UniProt -> Unified Protein Sequences Schema Mapping
  ============================================================
  Source: ./schema.json
  Target: ../schema.json (Unified Protein Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────┬─────────────────┬──────────────┐
  │ Source Field            │ Target Field    │ Transform    │
  ├─────────────────────────┼─────────────────┼──────────────┤
  │ primaryAccession        │ accession       │ Direct       │
  │ sequence.value          │ sequence        │ Direct       │
  │ sequence.length         │ length          │ Direct       │
  │ organism.scientificName │ organism        │ Direct       │
  │ organism.lineage        │ taxonomy        │ Direct       │
  │ genes[0].geneName       │ gene            │ Extract      │
  │ features                │ features        │ Restructure  │
  │ (computed)              │ _source         │ Metadata     │
  └─────────────────────────┴─────────────────┴──────────────┘

  NOTES:
  - UniProt is the definitive protein resource
  - Rich annotation maps to unified schema
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    exclude-result-prefixes="xs fn map array">

  <xsl:output method="json" indent="yes"/>

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- Primary accession -->
      <fn:string key="accession">
        <xsl:value-of select="fn:string[@key='primaryAccession']"/>
      </fn:string>

      <!-- Entry name -->
      <fn:string key="entry_name">
        <xsl:value-of select="fn:string[@key='uniProtkbId']"/>
      </fn:string>

      <!-- Reviewed status -->
      <fn:boolean key="reviewed">
        <xsl:value-of select="contains(fn:string[@key='entryType'], 'reviewed')"/>
      </fn:boolean>

      <!-- Sequence -->
      <fn:string key="sequence">
        <xsl:value-of select="fn:map[@key='sequence']/fn:string[@key='value']"/>
      </fn:string>

      <!-- Length -->
      <fn:number key="length">
        <xsl:value-of select="fn:map[@key='sequence']/fn:number[@key='length']"/>
      </fn:number>

      <!-- Molecular weight -->
      <xsl:if test="fn:map[@key='sequence']/fn:number[@key='molWeight']">
        <fn:number key="molecular_weight">
          <xsl:value-of select="fn:map[@key='sequence']/fn:number[@key='molWeight']"/>
        </fn:number>
      </xsl:if>

      <!-- Organism -->
      <fn:string key="organism">
        <xsl:value-of select="fn:map[@key='organism']/fn:string[@key='scientificName']"/>
      </fn:string>

      <!-- Taxonomy -->
      <xsl:if test="fn:map[@key='organism']/fn:array[@key='lineage']">
        <fn:array key="taxonomy">
          <xsl:for-each select="fn:map[@key='organism']/fn:array[@key='lineage']/fn:string">
            <fn:string><xsl:value-of select="."/></fn:string>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Gene name -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='genes']/fn:map[1]/fn:map[@key='geneName']/fn:string[@key='value']">
          <fn:string key="gene">
            <xsl:value-of select="fn:array[@key='genes']/fn:map[1]/fn:map[@key='geneName']/fn:string[@key='value']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Definition -->
      <xsl:if test="fn:map[@key='proteinDescription']/fn:map[@key='recommendedName']/fn:map[@key='fullName']/fn:string[@key='value']">
        <fn:string key="definition">
          <xsl:value-of select="fn:map[@key='proteinDescription']/fn:map[@key='recommendedName']/fn:map[@key='fullName']/fn:string[@key='value']"/>
        </fn:string>
      </xsl:if>

      <!-- Protein existence -->
      <fn:string key="protein_existence">
        <xsl:value-of select="fn:string[@key='proteinExistence']"/>
      </fn:string>

      <!-- Features -->
      <xsl:if test="fn:array[@key='features']">
        <fn:array key="features">
          <xsl:for-each select="fn:array[@key='features']/fn:map">
            <fn:map>
              <fn:string key="type">
                <xsl:value-of select="fn:string[@key='type']"/>
              </fn:string>
              <fn:map key="location">
                <fn:number key="start">
                  <xsl:value-of select="fn:map[@key='location']/fn:map[@key='start']/fn:number[@key='value']"/>
                </fn:number>
                <fn:number key="end">
                  <xsl:value-of select="fn:map[@key='location']/fn:map[@key='end']/fn:number[@key='value']"/>
                </fn:number>
              </fn:map>
              <xsl:if test="fn:string[@key='description']">
                <fn:string key="description">
                  <xsl:value-of select="fn:string[@key='description']"/>
                </fn:string>
              </xsl:if>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Cross-references -->
      <xsl:if test="fn:array[@key='uniProtKBCrossReferences']">
        <fn:array key="dbxrefs">
          <xsl:for-each select="fn:array[@key='uniProtKBCrossReferences']/fn:map">
            <fn:map>
              <fn:string key="database">
                <xsl:value-of select="fn:string[@key='database']"/>
              </fn:string>
              <fn:string key="id">
                <xsl:value-of select="fn:string[@key='id']"/>
              </fn:string>
            </fn:map>
          </xsl:for-each>
        </fn:array>
      </xsl:if>

      <!-- Source metadata -->
      <fn:map key="_source">
        <fn:string key="database">UniProt</fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='primaryAccession']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>

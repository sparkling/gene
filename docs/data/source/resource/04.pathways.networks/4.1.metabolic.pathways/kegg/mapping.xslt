<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  KEGG Pathway Database -> Unified Metabolic Pathways Schema Mapping
  ============================================================
  Source: ./schema.json (KEGG KGML/REST API data)
  Target: ../schema.json (Unified 4.1 Metabolic Pathways Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ pathway_id                  │ pathway_id              │ Direct       │
  │ name                        │ pathway_name            │ Direct       │
  │ org                         │ organism                │ Direct       │
  │ entry                       │ gene_entries            │ Restructure  │
  │ compound                    │ compound_entries        │ Restructure  │
  │ reaction                    │ reactions               │ Restructure  │
  │ link                        │ cross_references        │ Restructure  │
  │ entry_type                  │ kegg_entry_type         │ Direct       │
  │ relation_type               │ kegg_relation_type      │ Direct       │
  │ relation_subtype            │ kegg_relation_subtype   │ Direct       │
  │ graphics                    │ kegg_graphics           │ Direct       │
  │ reaction_type               │ kegg_reaction_type      │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - Missing optional fields -> null

  NOTES:
  - KEGG uses organism prefix (e.g., "hsa") in pathway IDs
  - Entry types follow KGML specification
  - Graphics coordinates use KGML positioning system
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== DIRECT FIELD MAPPINGS ========== -->

      <!-- Primary identifier (e.g., hsa00010) -->
      <fn:string key="pathway_id">
        <xsl:value-of select="fn:string[@key='pathway_id']"/>
      </fn:string>

      <!-- Pathway name -->
      <fn:string key="pathway_name">
        <xsl:value-of select="fn:string[@key='name']"/>
      </fn:string>

      <!-- Organism code (e.g., hsa, mmu) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='org'] and fn:string[@key='org'] != ''">
          <fn:string key="organism">
            <xsl:value-of select="fn:string[@key='org']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ARRAY TRANSFORMATIONS ========== -->

      <!-- Gene entries -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='entry']">
          <fn:array key="gene_entries">
            <xsl:for-each select="fn:array[@key='entry']/fn:map">
              <fn:map>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='id']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='name']">
                  <fn:string key="symbol">
                    <xsl:value-of select="fn:string[@key='name']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='type']">
                  <fn:string key="type">
                    <xsl:value-of select="fn:string[@key='type']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_entries"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Compound entries -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='compound']">
          <fn:array key="compound_entries">
            <xsl:for-each select="fn:array[@key='compound']/fn:map">
              <fn:map>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='id']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='name']">
                  <fn:string key="name">
                    <xsl:value-of select="fn:string[@key='name']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="compound_entries"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Reactions -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='reaction']">
          <fn:array key="reactions">
            <xsl:for-each select="fn:array[@key='reaction']/fn:map">
              <fn:map>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='id']"/>
                </fn:string>
                <xsl:if test="fn:array[@key='substrate']">
                  <fn:array key="substrates">
                    <xsl:for-each select="fn:array[@key='substrate']/fn:string">
                      <fn:string><xsl:value-of select="."/></fn:string>
                    </xsl:for-each>
                  </fn:array>
                </xsl:if>
                <xsl:if test="fn:array[@key='product']">
                  <fn:array key="products">
                    <xsl:for-each select="fn:array[@key='product']/fn:string">
                      <fn:string><xsl:value-of select="."/></fn:string>
                    </xsl:for-each>
                  </fn:array>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="reactions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Cross-references -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='link']">
          <fn:array key="cross_references">
            <xsl:for-each select="fn:array[@key='link']/fn:map">
              <fn:map>
                <fn:string key="database">
                  <xsl:value-of select="fn:string[@key='db']"/>
                </fn:string>
                <fn:string key="id">
                  <xsl:value-of select="fn:string[@key='id']"/>
                </fn:string>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cross_references"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== KEGG-SPECIFIC FIELDS ========== -->

      <!-- Entry type -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='entry_type'] and fn:string[@key='entry_type'] != ''">
          <fn:string key="kegg_entry_type">
            <xsl:value-of select="fn:string[@key='entry_type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="kegg_entry_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Relation type -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='relation_type'] and fn:string[@key='relation_type'] != ''">
          <fn:string key="kegg_relation_type">
            <xsl:value-of select="fn:string[@key='relation_type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="kegg_relation_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Relation subtype (array) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='relation_subtype']">
          <fn:array key="kegg_relation_subtype">
            <xsl:for-each select="fn:array[@key='relation_subtype']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="kegg_relation_subtype"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Graphics object -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='graphics']">
          <fn:map key="kegg_graphics">
            <xsl:copy-of select="fn:map[@key='graphics']/*"/>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="kegg_graphics"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Reaction type (reversible/irreversible) -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='reaction_type'] and fn:string[@key='reaction_type'] != ''">
          <fn:string key="kegg_reaction_type">
            <xsl:value-of select="fn:string[@key='reaction_type']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="kegg_reaction_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">KEGG</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='pathway_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Reusable Functions
       ============================================================ -->

  <!-- Convert empty/placeholder to empty string (for null handling) -->
  <xsl:function name="local:nullify" as="xs:string">
    <xsl:param name="value" as="xs:string"/>
    <xsl:choose>
      <xsl:when test="$value = ('-', '.', '', 'NA', 'N/A', 'null')"></xsl:when>
      <xsl:otherwise><xsl:value-of select="$value"/></xsl:otherwise>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>

<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  DrugBank -> Pharmaceuticals Unified Schema Mapping
  ============================================================
  Source: ./schema.md (DrugBank Comprehensive Drug Database)
  Target: ../schema.json (2.2 Pharmaceuticals Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ drugbank_id                 │ drug_id             │ Direct       │
  │ drugbank_id                 │ drugbank_id         │ Direct       │
  │ name                        │ name                │ Direct       │
  │ smiles                      │ canonical_smiles    │ Direct       │
  │ inchikey                    │ inchi_key           │ Direct       │
  │ average_mass                │ molecular_weight    │ Direct       │
  │ chemical_formula            │ molecular_formula   │ Direct       │
  │ drug_type                   │ drug_type           │ Direct       │
  │ status                      │ approval_status     │ Direct       │
  │ cas_number                  │ cas_number          │ Direct       │
  │ mechanism_of_action         │ mechanism_of_action │ Direct       │
  │ drug_interactions           │ drug_interactions   │ Object array │
  │ half_life                   │ half_life           │ Direct       │
  │ protein_binding             │ protein_binding     │ Direct       │
  │ targets                     │ targets             │ Object array │
  │ pathways                    │ pathways            │ Object array │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null
  - status: maps "approved", "investigational", "experimental", etc.

  NOTES:
  - DrugBank uses DB prefix for IDs (e.g., DB00945)
  - Comprehensive drug data including PK/PD, interactions, pathways
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

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== IDENTIFIERS ========== -->

      <fn:string key="drug_id">
        <xsl:value-of select="fn:string[@key='drugbank_id']"/>
      </fn:string>

      <fn:string key="drugbank_id">
        <xsl:value-of select="fn:string[@key='drugbank_id']"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='cas_number'])">
          <fn:null key="cas_number"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="cas_number">
            <xsl:value-of select="fn:string[@key='cas_number']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='smiles'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchikey'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='inchikey']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MOLECULAR PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='chemical_formula'])">
          <fn:null key="molecular_formula"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_formula">
            <xsl:value-of select="fn:string[@key='chemical_formula']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="not(fn:number[@key='average_mass']) or fn:number[@key='average_mass'] &lt;= 0">
          <fn:null key="molecular_weight"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='average_mass']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG CLASSIFICATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='drug_type'])">
          <fn:null key="drug_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="drug_type">
            <xsl:value-of select="fn:string[@key='drug_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='status'])">
          <fn:null key="approval_status"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="approval_status">
            <xsl:value-of select="fn:string[@key='status']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PHARMACOLOGY ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='mechanism_of_action'])">
          <fn:null key="mechanism_of_action"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="mechanism_of_action">
            <xsl:value-of select="fn:string[@key='mechanism_of_action']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='half_life'])">
          <fn:null key="half_life"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="half_life">
            <xsl:value-of select="fn:string[@key='half_life']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='protein_binding'])">
          <fn:null key="protein_binding"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="protein_binding">
            <xsl:value-of select="fn:string[@key='protein_binding']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG INTERACTIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='drug_interactions']/fn:map">
          <fn:array key="drug_interactions">
            <xsl:for-each select="fn:array[@key='drug_interactions']/fn:map">
              <fn:map>
                <fn:string key="drug">
                  <xsl:value-of select="fn:string[@key='drug']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='description']">
                  <fn:string key="description">
                    <xsl:value-of select="fn:string[@key='description']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='severity']">
                  <fn:string key="severity">
                    <xsl:value-of select="fn:string[@key='severity']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="drug_interactions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== TARGETS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='targets']/fn:map">
          <fn:array key="targets">
            <xsl:for-each select="fn:array[@key='targets']/fn:map">
              <fn:map>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='name']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='uniprot_id']">
                  <fn:string key="uniprot_id">
                    <xsl:value-of select="fn:string[@key='uniprot_id']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='action']">
                  <fn:string key="action">
                    <xsl:value-of select="fn:string[@key='action']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="targets"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PATHWAYS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='pathways']/fn:map">
          <fn:array key="pathways">
            <xsl:for-each select="fn:array[@key='pathways']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='pathway_id']">
                  <fn:string key="pathway_id">
                    <xsl:value-of select="fn:string[@key='pathway_id']"/>
                  </fn:string>
                </xsl:if>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='name']"/>
                </fn:string>
                <xsl:if test="fn:string[@key='source']">
                  <fn:string key="source">
                    <xsl:value-of select="fn:string[@key='source']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pathways"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">DrugBank</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='drugbank_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>

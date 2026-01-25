<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  HPO Annotations -> Phenotype Databases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (HPO Disease-Phenotype Annotations)
  Target: ../schema.json (3.2 Phenotype Databases Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ DatabaseID                  │ disease_id              │ Direct       │
  │ DiseaseName                 │ disease_name            │ Direct       │
  │ HPO_ID                      │ hpo_ids                 │ Array        │
  │ Frequency                   │ frequency               │ Direct       │
  │ Gene_ID                     │ gene_id                 │ Direct       │
  │ Gene_Symbol                 │ gene_symbol             │ Direct       │
  │ Qualifier                   │ qualifier               │ Direct       │
  │ Reference                   │ reference               │ Direct       │
  │ Evidence                    │ evidence                │ Enum         │
  │ Onset                       │ onset                   │ Direct       │
  │ Sex                         │ sex                     │ Enum         │
  │ Modifier                    │ modifier                │ Array        │
  │ Aspect                      │ aspect                  │ Enum         │
  │ Biocuration                 │ biocuration             │ Direct       │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - HPO annotations link diseases to phenotypes
  - Evidence codes: IEA, TAS, PCS
  - Aspect: P (Phenotype), I (Inheritance), C (Course), M (Modifier)
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
      <!-- ========== DISEASE IDENTIFIERS ========== -->

      <fn:string key="disease_id">
        <xsl:value-of select="fn:string[@key='DatabaseID']"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='DiseaseName'])">
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='DiseaseName']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== HPO PHENOTYPES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='HPO_ID']/fn:string">
          <fn:array key="hpo_ids">
            <xsl:for-each select="fn:array[@key='HPO_ID']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='HPO_ID'] and not(local:is-null(fn:string[@key='HPO_ID']))">
          <fn:array key="hpo_ids">
            <fn:string><xsl:value-of select="fn:string[@key='HPO_ID']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="hpo_ids"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Frequency'])">
          <fn:null key="frequency"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="frequency">
            <xsl:value-of select="fn:string[@key='Frequency']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE INFORMATION ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='Gene_ID']">
          <fn:number key="gene_id">
            <xsl:value-of select="fn:number[@key='Gene_ID']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="gene_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Gene_Symbol'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='Gene_Symbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ANNOTATION ATTRIBUTES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Qualifier'])">
          <fn:null key="qualifier"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="qualifier">
            <xsl:value-of select="fn:string[@key='Qualifier']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Reference'])">
          <fn:null key="reference"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="reference">
            <xsl:value-of select="fn:string[@key='Reference']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Evidence'])">
          <fn:null key="evidence"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="evidence">
            <xsl:value-of select="fn:string[@key='Evidence']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Onset'])">
          <fn:null key="onset"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="onset">
            <xsl:value-of select="fn:string[@key='Onset']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Sex'])">
          <fn:null key="sex"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="sex">
            <xsl:value-of select="fn:string[@key='Sex']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:array[@key='Modifier']/fn:string">
          <fn:array key="modifier">
            <xsl:for-each select="fn:array[@key='Modifier']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="modifier"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Aspect'])">
          <fn:null key="aspect"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="aspect">
            <xsl:value-of select="fn:string[@key='Aspect']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Biocuration'])">
          <fn:null key="biocuration"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="biocuration">
            <xsl:value-of select="fn:string[@key='Biocuration']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">HPO Annotations</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='DatabaseID']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>

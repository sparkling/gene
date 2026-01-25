<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Pathway Commons -> Unified Signaling Pathways Schema Mapping
  ============================================================
  Source: ./schema.json (Pathway Commons BioPAX data)
  Target: ../schema.json (Unified 4.2 Signaling Pathways Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ uri                         │ uri                     │ Direct       │
  │ displayName                 │ display_name            │ Direct       │
  │ organism                    │ organism                │ Object       │
  │ dataSource                  │ data_source             │ Direct       │
  │ pathwayComponent            │ pathway_components      │ Array        │
  │ pathwayOrder                │ pathway_order           │ Array        │
  │ left                        │ left                    │ Array        │
  │ right                       │ right                   │ Array        │
  │ conversionDirection         │ conversion_direction    │ Direct       │
  │ eCNumber                    │ ec_numbers              │ Array        │
  │ controlType                 │ control_type            │ Direct       │
  │ controller                  │ controller              │ Array        │
  │ controlled                  │ controlled              │ Direct       │
  │ entityReference             │ entity_reference        │ Object       │
  │ feature                     │ features                │ Array        │
  │ cellularLocation            │ cellular_location       │ Object       │
  │ component                   │ components              │ Array        │
  │ componentStoichiometry      │ component_stoichiometry │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NULL HANDLING:
  - Empty strings "" -> null
  - Missing optional fields -> null

  NOTES:
  - Pathway Commons uses BioPAX URIs as identifiers
  - Aggregates data from Reactome, KEGG, WikiPathways, PID, PhosphoSitePlus, etc.
  - conversionDirection uses LEFT-TO-RIGHT, RIGHT-TO-LEFT, REVERSIBLE
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
      <!-- ========== REQUIRED FIELDS ========== -->

      <!-- BioPAX URI -->
      <fn:string key="uri">
        <xsl:value-of select="fn:string[@key='uri']"/>
      </fn:string>

      <!-- Display name -->
      <fn:string key="display_name">
        <xsl:value-of select="fn:string[@key='displayName']"/>
      </fn:string>

      <!-- ========== ORGANISM ========== -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='organism']">
          <fn:map key="organism">
            <xsl:if test="fn:map[@key='organism']/fn:string[@key='displayName']">
              <fn:string key="display_name">
                <xsl:value-of select="fn:map[@key='organism']/fn:string[@key='displayName']"/>
              </fn:string>
            </xsl:if>
            <xsl:if test="fn:map[@key='organism']/fn:number[@key='taxId']">
              <fn:number key="tax_id">
                <xsl:value-of select="fn:map[@key='organism']/fn:number[@key='taxId']"/>
              </fn:number>
            </xsl:if>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Data source -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='dataSource'] and fn:string[@key='dataSource'] != ''">
          <fn:string key="data_source">
            <xsl:value-of select="fn:string[@key='dataSource']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="data_source"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PATHWAY STRUCTURE ========== -->

      <!-- Pathway components -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='pathwayComponent']">
          <fn:array key="pathway_components">
            <xsl:for-each select="fn:array[@key='pathwayComponent']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pathway_components"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Pathway order -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='pathwayOrder']">
          <fn:array key="pathway_order">
            <xsl:for-each select="fn:array[@key='pathwayOrder']/fn:map">
              <fn:map>
                <xsl:copy-of select="./*"/>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="pathway_order"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== REACTION FIELDS ========== -->

      <!-- Left (substrates) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='left']">
          <fn:array key="left">
            <xsl:for-each select="fn:array[@key='left']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="left"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Right (products) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='right']">
          <fn:array key="right">
            <xsl:for-each select="fn:array[@key='right']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="right"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Conversion direction -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='conversionDirection'] and fn:string[@key='conversionDirection'] != ''">
          <fn:string key="conversion_direction">
            <xsl:value-of select="fn:string[@key='conversionDirection']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="conversion_direction"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- EC numbers -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='eCNumber']">
          <fn:array key="ec_numbers">
            <xsl:for-each select="fn:array[@key='eCNumber']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ec_numbers"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CONTROL FIELDS ========== -->

      <!-- Control type -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='controlType'] and fn:string[@key='controlType'] != ''">
          <fn:string key="control_type">
            <xsl:value-of select="fn:string[@key='controlType']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="control_type"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Controller -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='controller']">
          <fn:array key="controller">
            <xsl:for-each select="fn:array[@key='controller']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="controller"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Controlled process -->
      <xsl:choose>
        <xsl:when test="fn:string[@key='controlled'] and fn:string[@key='controlled'] != ''">
          <fn:string key="controlled">
            <xsl:value-of select="fn:string[@key='controlled']"/>
          </fn:string>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="controlled"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ENTITY FIELDS ========== -->

      <!-- Entity reference -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='entityReference']">
          <fn:map key="entity_reference">
            <xsl:if test="fn:map[@key='entityReference']/fn:string[@key='db']">
              <fn:string key="db">
                <xsl:value-of select="fn:map[@key='entityReference']/fn:string[@key='db']"/>
              </fn:string>
            </xsl:if>
            <xsl:if test="fn:map[@key='entityReference']/fn:string[@key='id']">
              <fn:string key="id">
                <xsl:value-of select="fn:map[@key='entityReference']/fn:string[@key='id']"/>
              </fn:string>
            </xsl:if>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="entity_reference"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Features (PTMs, domains) -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='feature']">
          <fn:array key="features">
            <xsl:for-each select="fn:array[@key='feature']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='type']">
                  <fn:string key="type">
                    <xsl:value-of select="fn:string[@key='type']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='position']">
                  <fn:string key="position">
                    <xsl:value-of select="fn:string[@key='position']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="features"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Cellular location -->
      <xsl:choose>
        <xsl:when test="fn:map[@key='cellularLocation']">
          <fn:map key="cellular_location">
            <xsl:if test="fn:map[@key='cellularLocation']/fn:string[@key='term']">
              <fn:string key="term">
                <xsl:value-of select="fn:map[@key='cellularLocation']/fn:string[@key='term']"/>
              </fn:string>
            </xsl:if>
            <xsl:if test="fn:map[@key='cellularLocation']/fn:string[@key='go']">
              <fn:string key="go">
                <xsl:value-of select="fn:map[@key='cellularLocation']/fn:string[@key='go']"/>
              </fn:string>
            </xsl:if>
          </fn:map>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="cellular_location"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== COMPLEX FIELDS ========== -->

      <!-- Components -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='component']">
          <fn:array key="components">
            <xsl:for-each select="fn:array[@key='component']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="components"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Component stoichiometry -->
      <xsl:choose>
        <xsl:when test="fn:array[@key='componentStoichiometry']">
          <fn:array key="component_stoichiometry">
            <xsl:for-each select="fn:array[@key='componentStoichiometry']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='entity']">
                  <fn:string key="entity">
                    <xsl:value-of select="fn:string[@key='entity']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:number[@key='coefficient']">
                  <fn:number key="coefficient">
                    <xsl:value-of select="fn:number[@key='coefficient']"/>
                  </fn:number>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="component_stoichiometry"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Pathway Commons</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='uri']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

</xsl:stylesheet>

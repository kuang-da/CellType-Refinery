// @ts-check
import {themes as prismThemes} from 'prism-react-renderer';

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'CellType-Refinery',
  tagline: 'Cell-type annotation for spatial proteomics',
  favicon: 'img/favicon.ico',

  url: 'https://kuang-da.github.io',
  baseUrl: '/CellType-Refinery/',
  trailingSlash: false,

  organizationName: 'kuang-da',
  projectName: 'CellType-Refinery',

  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          sidebarPath: './sidebars.js',
          editUrl: 'https://github.com/kuang-da/CellType-Refinery/tree/main/celltype-refinery-docs/',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      }),
    ],
  ],

  themes: ['@docusaurus/theme-mermaid'],

  markdown: {
    mermaid: true,
  },

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      image: 'img/social-card.png',
      navbar: {
        title: 'CellType-Refinery',
        logo: {
          alt: 'CellType-Refinery Logo',
          src: 'img/logo.svg',
        },
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'docsSidebar',
            position: 'left',
            label: 'Documentation',
          },
          {
            href: 'https://github.com/kuang-da/CellType-Refinery',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Documentation',
            items: [
              { label: 'Getting Started', to: '/docs/getting-started/installation' },
              { label: 'Core Workflows', to: '/docs/core-workflows/workflow-overview' },
              { label: 'CLI Reference', to: '/docs/cli/overview' },
            ],
          },
          {
            title: 'More',
            items: [
              { label: 'GitHub', href: 'https://github.com/kuang-da/CellType-Refinery' },
              { label: 'Issues', href: 'https://github.com/kuang-da/CellType-Refinery/issues' },
            ],
          },
        ],
        copyright: `Copyright ${new Date().getFullYear()} Penn Biomedical Image Analysis Lab. Built with Docusaurus.`,
      },
      prism: {
        theme: prismThemes.github,
        darkTheme: prismThemes.dracula,
        additionalLanguages: ['python', 'bash', 'yaml', 'json'],
      },
      mermaid: {
        theme: {light: 'neutral', dark: 'dark'},
      },
    }),
};

export default config;
